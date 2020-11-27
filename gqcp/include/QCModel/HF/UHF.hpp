// This file is part of GQCG-GQCP.
//
// Copyright (C) 2017-2020  the GQCG developers
//
// GQCG-GQCP is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GQCG-GQCP is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-GQCP.  If not, see <http://www.gnu.org/licenses/>.

#pragma once


#include "Basis/SpinorBasis/SpinResolvedOrbitalSpace.hpp"
#include "Basis/Transformations/UTransformation.hpp"
#include "Basis/Transformations/UTransformationComponent.hpp"
#include "DensityMatrix/SpinResolved1DM.hpp"
#include "DensityMatrix/SpinResolved2DM.hpp"
#include "Mathematical/Representation/Matrix.hpp"
#include "Operator/SecondQuantized/RSQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/USQOneElectronOperator.hpp"
#include "QCModel/HF/RHF.hpp"
#include "QuantumChemical/Spin.hpp"


namespace GQCP {
namespace QCModel {


/**
 *  The unrestricted Hartree-Fock wave function model.
 * 
 *  @tparam _Scalar         the scalar representation of the coefficients in the coefficient matrix
 */
template <typename _Scalar>
class UHF {
public:
    using Scalar = _Scalar;


private:
    size_t N_alpha;  // the number of alpha electrons
    size_t N_beta;   // the number of beta electrons

    VectorX<double> orbital_energies_alpha;  // sorted by ascending energy
    VectorX<double> orbital_energies_beta;   // sorted by ascending energy

    // The transformation between the UHF MOs and the atomic spin-orbitals.
    UTransformation<Scalar> C;


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  A basic constructor that sets all the properties.
     * 
     *  @param N_alpha                                  the number of alpha electrons
     *  @param N_beta                                   the number of beta electrons
     *  @param orbital_energies_alpha                   the orbital energies for the alpha-spin-orbitals, sorted by ascending energy
     *  @param orbital_energies_beta                    the orbital energies for the beta-spin-orbitals, sorted by ascending energy
     *  @param C                                        The transformation between the UHF MOs and the atomic spin-orbitals.
     */
    UHF(const size_t N_alpha, const size_t N_beta, const VectorX<double>& orbital_energies_alpha, const VectorX<double>& orbital_energies_beta, const UTransformation<Scalar>& C) :
        N_alpha {N_alpha},
        N_beta {N_beta},
        orbital_energies_alpha {orbital_energies_alpha},
        orbital_energies_beta {orbital_energies_beta},
        C {C} {

        // Check for valid arguments.
        const auto K_alpha = C.alpha().numberOfOrbitals();  // number of alpha spatial orbitals
        const auto K_beta = C.beta().numberOfOrbitals();    // number of beta spatial orbitals

        if (N_alpha > K_alpha) {
            throw std::invalid_argument("UHF(const size_t, const size_t, const VectorX<double>&, const VectorX<double>&, const UTransformation<Scalar>&): The number of given alpha electrons cannot be larger than the number of alpha spatial orbitals.");
        }

        if (N_beta > K_beta) {
            throw std::invalid_argument("UHF(const size_t, const size_t, const VectorX<double>&, const VectorX<double>&, const UTransformation<Scalar>&): The number of given beta electrons cannot be larger than the number of beta spatial orbitals.");
        }

        if (K_alpha != orbital_energies_alpha.size()) {
            throw std::invalid_argument("UHF(const size_t, const size_t, const VectorX<double>&, const VectorX<double>&, const UTransformation<Scalar>&): The number of given alpha-spin-orbital energies does not match the number of alpha spin-orbitals.");
        }

        if (K_beta != orbital_energies_beta.size()) {
            throw std::invalid_argument("UHF(const size_t, const size_t, const VectorX<double>&, const VectorX<double>&, const UTransformation<Scalar>&): The number of given beta-spin-orbital energies does not match the beta of beta spin-orbitals.");
        }
    }


    /**
     *  A constructor that initializes the orbital energies to zeros.
     * 
     *  @param N_alpha                                  the number of alpha electrons
     *  @param N_beta                                   the number of beta electrons
     *  @param orbital_energies_alpha                   the orbital energies for the alpha-spin-orbitals, sorted by ascending energy
     *  @param orbital_energies_beta                    the orbital energies for the beta-spin-orbitals, sorted by ascending energy
     *  @param C                                        The transformation between the UHF MOs and the atomic spin-orbitals.
     */
    UHF(const size_t N_alpha, const size_t N_beta, const UTransformation<Scalar>& C) :
        UHF(N_alpha, N_beta,
            GQCP::VectorX<double>::Zero(C.component(Spin::alpha).numberOfOrbitals()),
            GQCP::VectorX<double>::Zero(C.component(Spin::beta).numberOfOrbitals()),
            C) {
    }


    /**
     *  Convert an RHF wave function model to an UHF wave function model.
     * 
     *  @param rhf_model            an RHF wave function model
     */
    UHF(const GQCP::QCModel::RHF<Scalar>& rhf_model) :
        UHF(rhf_model.numberOfElectrons(Spin::alpha), rhf_model.numberOfElectrons(Spin::beta),
            rhf_model.orbitalEnergies(), rhf_model.orbitalEnergies(),
            GQCP::UTransformation<Scalar>::FromRestricted(rhf_model.expansion())) {}


    /*
     *  PUBLIC STATIC METHODS
     */

    /**
     *  @param P_sigma              the spin-sigma density matrix in a scalar basis
     *  @param H_core_sigma         the spin-sigma core Hamiltonian expressed in the same scalar basis
     *  @param F_sigma              the spin-sigma Fock matrix in the same scalar basis
     *
     *  @return the UHF electronic energy for the sigma electrons
     */
    static double calculateElectronicEnergy(const SpinResolved1DMComponent<Scalar>& P_sigma, const ScalarRSQOneElectronOperator<Scalar>& H_core_sigma, const ScalarRSQOneElectronOperator<Scalar>& F_sigma) {

        // First, calculate the sum of H_core and F (this saves a contraction).
        const auto Z_sigma = H_core_sigma + F_sigma;

        // Convert the matrix Z to a GQCP::Tensor<double, 2> Z_tensor.
        // Einsum is only implemented for a tensor + a matrix, not for 2 matrices.
        Eigen::TensorMap<Eigen::Tensor<const Scalar, 2>> Z_sigma_t {Z_sigma.parameters().data(), Z_sigma.parameters().rows(), Z_sigma.parameters().cols()};
        Tensor<Scalar, 2> Z_sigma_tensor = Tensor<Scalar, 2>(Z_sigma_t);

        // To calculate the electronic energy, we must perform a double contraction (with prefactor 0.5).
        //      0.5 P_sigma(mu nu) P_sigma(mu nu)
        Tensor<Scalar, 0> contraction = 0.5 * Z_sigma_tensor.template einsum<2>("ij,ji->", P_sigma);

        // As the double contraction of two rank-2 tensors is a scalar (a tensor of rank 0), we should access the value as (0).
        return contraction(0);
    }


    /**
     *  @param F_sigma                  the sigma-spin Fock matrix expressed in the AO basis
     *  @param D_sigma                  the sigma-spin density matrix in the AO basis
     *  @param S                        the overlap matrix of the AO basis
     * 
     *  @return the sigma-spin error matrix
     */
    static SquareMatrix<Scalar> calculateError(const SquareMatrix<Scalar>& F_sigma, const SpinResolved1DMComponent<Scalar>& D_sigma, const SquareMatrix<Scalar>& S) {
        return QCModel::RHF<Scalar>::calculateError(F_sigma, D_sigma, S);
    }


    /**
     *  Calculate and return the (spin-resolved) UHF 1-DM in an orthonormal spin-orbital basis.
     * 
     *  @param K_a          The number of spatial orbitals for the alpha spin component.
     *  @param K_b          The number of spatial orbitals for the beta spin component.
     *  @param N_a          The number of alpha electrons, i.e. the number of occupied alpha spin-orbitals.
     *  @param N_b          The number of beta electrons, i.e. the number of occupied beta spin-orbitals.
     *
     *  @return The UHF 1-DM in an orthonormal spin-orbital basis.
     */
    static SpinResolved1DM<Scalar> calculateOrthonormalBasis1DM(const size_t K_a, const size_t K_b, const size_t N_a, const size_t N_b) {

        // The (alpha) 1-DM for UHF looks like (for K_alpha=5, N_alpha=3)
        //    1  0  0  0  0
        //    0  1  0  0  0
        //    0  0  1  0  0
        //    0  0  0  0  0
        //    0  0  0  0  0

        SpinResolved1DMComponent<Scalar> D_MO_a = SpinResolved1DMComponent<Scalar>::Zero(K_a);
        D_MO_a.topLeftCorner(N_a, N_a) = SquareMatrix<Scalar>::Identity(N_a);

        SpinResolved1DMComponent<Scalar> D_MO_b = SpinResolved1DMComponent<Scalar>::Zero(K_b);
        D_MO_b.topLeftCorner(N_b, N_b) = SquareMatrix<Scalar>::Identity(N_b);

        return SpinResolved1DM<Scalar> {D_MO_a, D_MO_b};
    }


    /**
     *  Calculate and return the (spin-resolved) UHF 2-DM in an orthonormal spin-orbital basis.
     * 
     *  @param K            The number of spatial orbitals for the alpha and beta spin components.
     *  @param N_a          The number of alpha electrons, i.e. the number of occupied alpha spin-orbitals.
     *  @param N_b          The number of beta electrons, i.e. the number of occupied beta spin-orbitals.
     *
     *  @return The UHF 2-DM in an orthonormal spin-orbital basis.
     */
    static SpinResolved2DM<Scalar> calculateOrthonormalBasis2DM(const size_t K, const size_t N_a, const size_t N_b) {

        // Create the orbital space to determine the loops.
        const auto orbital_space = UHF<Scalar>::orbitalSpace(K, K, N_a, N_b);


        // Use KISS formulas to implement the spin components of the UHF 2-DM.
        PureSpinResolved2DMComponent<Scalar> d_aaaa = PureSpinResolved2DMComponent<Scalar>::Zero(K);
        for (const auto& i : orbital_space.alpha().indices(OccupationType::k_occupied)) {
            for (const auto& j : orbital_space.alpha().indices(OccupationType::k_occupied)) {
                for (const auto& k : orbital_space.alpha().indices(OccupationType::k_occupied)) {
                    for (const auto& l : orbital_space.alpha().indices(OccupationType::k_occupied)) {
                        if ((i == j) && (k == l)) {
                            d_aaaa(i, j, k, l) += 1.0;
                        }

                        if ((i == l) && (j == k)) {
                            d_aaaa(i, j, k, l) -= 1.0;
                        }
                    }
                }
            }
        }


        MixedSpinResolved2DMComponent<Scalar> d_aabb = MixedSpinResolved2DMComponent<Scalar>::Zero(K);
        for (const auto& i : orbital_space.alpha().indices(OccupationType::k_occupied)) {
            for (const auto& j : orbital_space.alpha().indices(OccupationType::k_occupied)) {
                for (const auto& k : orbital_space.beta().indices(OccupationType::k_occupied)) {
                    for (const auto& l : orbital_space.beta().indices(OccupationType::k_occupied)) {
                        if ((i == j) && (k == l)) {
                            d_aabb(i, j, k, l) += 1.0;
                        }
                    }
                }
            }
        }


        MixedSpinResolved2DMComponent<Scalar> d_bbaa = MixedSpinResolved2DMComponent<Scalar>::Zero(K);
        for (const auto& i : orbital_space.beta().indices(OccupationType::k_occupied)) {
            for (const auto& j : orbital_space.beta().indices(OccupationType::k_occupied)) {
                for (const auto& k : orbital_space.alpha().indices(OccupationType::k_occupied)) {
                    for (const auto& l : orbital_space.alpha().indices(OccupationType::k_occupied)) {
                        if ((i == j) && (k == l)) {
                            d_aabb(i, j, k, l) += 1.0;
                        }
                    }
                }
            }
        }


        PureSpinResolved2DMComponent<Scalar> d_bbbb = PureSpinResolved2DMComponent<Scalar>::Zero(K);
        for (const auto& i : orbital_space.beta().indices(OccupationType::k_occupied)) {
            for (const auto& j : orbital_space.beta().indices(OccupationType::k_occupied)) {
                for (const auto& k : orbital_space.beta().indices(OccupationType::k_occupied)) {
                    for (const auto& l : orbital_space.beta().indices(OccupationType::k_occupied)) {
                        if ((i == j) && (k == l)) {
                            d_aaaa(i, j, k, l) += 1.0;
                        }

                        if ((i == l) && (j == k)) {
                            d_aaaa(i, j, k, l) -= 1.0;
                        }
                    }
                }
            }
        }

        return SpinResolved2DM<Scalar> {d_aaaa, d_aabb, d_bbaa, d_bbbb};
    }


    /**
     *  Calculate the UHF 1-DM expressed in the underlying scalar basis.
     * 
     *  @param C            The transformation between the UHF MOs and the atomic spin-orbitals.
     *  @param N_a          The number of alpha electrons, i.e. the number of occupied alpha spin-orbitals.
     *  @param N_b          The number of beta electrons, i.e. the number of occupied beta spin-orbitals.
     *
     *  @return The UHF 1-DM expressed in the underlying scalar basis.
     */
    static SpinResolved1DM<Scalar> calculateScalarBasis1DM(const UTransformation<Scalar>& C, const size_t N_a, const size_t N_b) {

        // Calculate the 1-DM in the spin-orbital basis, and transform to the underlying scalar basis.
        const auto K_a = C.alpha().numberOfOrbitals();
        const auto K_b = C.beta().numberOfOrbitals();
        const auto D_orthonormal = UHF<Scalar>::calculateOrthonormalBasis1DM(K_a, K_b, N_a, N_b);

        return D_orthonormal.transformed(C.inverse());
    }


    /**
     *  Calculate the UHF direct (Coulomb) matrix for spin sigma.
     * 
     *  @param P_alpha              the UHF alpha density matrix expressed in the underlying scalar orbital basis
     *  @param P_beta               the UHF beta density matrix expressed in the (same) underlying scalar orbital basis
     *  @param sq_hamiltonian       the Hamiltonian expressed in the (same) underlying scalar orbital basis
     * 
     *  @return the UHF direct (Coulomb) matrix for spin sigma
     */
    static ScalarUSQOneElectronOperator<Scalar> calculateScalarBasisDirectMatrix(const SpinResolved1DM<Scalar>& P, const RSQHamiltonian<Scalar>& sq_hamiltonian) {

        // Get the two-electron parameters
        const auto& g = sq_hamiltonian.twoElectron().parameters();

        // Specify the contraction pairs for the direct contractions:
        //      (mu nu|rho lambda) P(rho lambda)
        const auto J_alpha = g.template einsum<2>("ijkl,kl->ij", P.alpha()).asMatrix();
        const auto J_beta = g.template einsum<2>("ijkl,kl->ij", P.beta()).asMatrix();

        // Calculate the total J tensor
        const auto J = J_alpha + J_beta;

        return ScalarUSQOneElectronOperator<Scalar> {J, J};
    }


    /**
     *  Calculate the UHF exchange matrix for spin sigma.
     * 
     *  @param P_sigma              the UHF sigma density matrix expressed in the underlying scalar orbital basis
     *  @param sq_hamiltonian       the Hamiltonian expressed in the (same) underlying scalar orbital basis
     * 
     *  @return the UHF direct (Coulomb) matrix for spin sigma
     */
    static ScalarUSQOneElectronOperator<Scalar> calculateScalarBasisExchangeMatrix(const SpinResolved1DM<Scalar>& P, const RSQHamiltonian<Scalar>& sq_hamiltonian) {

        // Get the two-electron parameters
        const auto& g = sq_hamiltonian.twoElectron().parameters();

        // Specify the contraction pairs for the exchange contraction:
        //      (mu rho|lambda nu) P(lambda rho)
        const auto K_alpha = g.template einsum<2>("ijkl,kj->il", P.alpha()).asMatrix();
        const auto K_beta = g.template einsum<2>("ijkl,kj->il", P.beta()).asMatrix();

        return ScalarUSQOneElectronOperator<Scalar> {K_alpha, K_beta};
    }


    /**
     *  Calculate the UHF Fock matrix F = H_core + G, in which G is a contraction of the density matrix and the two-electron integrals
     *
     *  @param P                    the RHF density matrix in a scalar basis
     *  @param sq_hamiltonian       the Hamiltonian expressed in the same scalar basis
     *
     *  @return the RHF Fock matrix expressed in the scalar basis
     */
    static ScalarUSQOneElectronOperator<Scalar> calculateScalarBasisFockMatrix(const SpinResolved1DM<Scalar>& P, const RSQHamiltonian<Scalar>& sq_hamiltonian) {

        // F_sigma = H_core + (J_alpha + J_beta) - K_sigma
        // H_core is always the same
        const auto& H_core = sq_hamiltonian.core().parameters();

        // Get the alpha and beta parameters of the coulomb and exchange matrices
        const auto J_a = UHF<Scalar>::calculateScalarBasisDirectMatrix(P, sq_hamiltonian).alpha().parameters();
        const auto J_b = UHF<Scalar>::calculateScalarBasisDirectMatrix(P, sq_hamiltonian).beta().parameters();

        const auto K_a = UHF<Scalar>::calculateScalarBasisExchangeMatrix(P, sq_hamiltonian).alpha().parameters();
        const auto K_b = UHF<Scalar>::calculateScalarBasisExchangeMatrix(P, sq_hamiltonian).beta().parameters();


        // Generate the alpha and beta Fock matrix and put the in a USQOneElectronOperator
        const auto F_a = H_core + J_a - K_a;
        const auto F_b = H_core + J_b - K_b;

        return ScalarUSQOneElectronOperator<Scalar> {F_a, F_b};
    }


    /**
     *  @param K_a            The number of alpha spatial orbitals.
     *  @param K_b            The number of beta spatial orbitals.
     *  @param N_a            The number of alpha electrons.
     *  @param N_b            The number of beta electrons.
     * 
     *  @return The implicit (i.e. with ascending and contiguous orbital indices) occupied-virtual orbital space for both the alpha and beta components that corresponds to these UHF model parameters.
     */
    static SpinResolvedOrbitalSpace orbitalSpace(const size_t K_a, const size_t K_b, const size_t N_a, const size_t N_b) {

        const auto alpha_orbital_space = OrbitalSpace::Implicit({{OccupationType::k_occupied, N_a}, {OccupationType::k_virtual, K_a - N_a}});
        const auto beta_orbital_space = OrbitalSpace::Implicit({{OccupationType::k_occupied, N_b}, {OccupationType::k_virtual, K_b - N_b}});

        return SpinResolvedOrbitalSpace {alpha_orbital_space, beta_orbital_space};
    }


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return The (spin-resolved) UHF 1-DM expressed in an orthonormal spin-orbital basis for these UHF model parameters.
     */
    SpinResolved1DM<Scalar> calculateOrthonormalBasis1DM() const {

        const auto K_a = this->numberOfSpinOrbitals(Spin::alpha);
        const auto K_b = this->numberOfSpinOrbitals(Spin::beta);
        const auto N_a = this->numberOfElectrons(Spin::alpha);
        const auto N_b = this->numberOfElectrons(Spin::beta);

        return UHF<Scalar>::calculateOrthonormalBasis1DM(K_a, K_b, N_a, N_b);
    }


    /**
     *  @return The (spin-resolved) UHF 2-DM expressed in an orthonormal spin-orbital basis for these UHF model parameters.
     * 
     *  @note We assume that the total number of alpha and beta orbitals is the same.
     */
    SpinResolved2DM<Scalar> calculateOrthonormalBasis2DM() const {

        const auto K = this->numberOfSpinOrbitals(Spin::alpha);  // Assume K_alpha and K_beta are equal.
        const auto N_a = this->numberOfElectrons(Spin::alpha);
        const auto N_b = this->numberOfElectrons(Spin::beta);

        return UHF<Scalar>::calculateOrthonormalBasis2DM(K, N_a, N_b);
    }


    /**
     *  @return The spin resolved UHF 1-DM expressed in the underlying scalar basis for these UHF model parameters.
     */
    SpinResolved1DM<Scalar> calculateScalarBasis1DM() const {

        const auto C_a = this->expansion().component(Spin::alpha);
        const auto C_b = this->expansion().component(Spin::beta);
        const UTransformation<Scalar> C {C_a, C_b};

        const auto N_a = this->numberOfElectrons(Spin::alpha);
        const auto N_b = this->numberOfElectrons(Spin::beta);

        return UHF<Scalar>::calculateScalarBasis1DM(C, N_a, N_b);
    }


    /**
     *  @return The transformation that expresses the UHF MOs in terms of the underlying AOs.
     */
    const UTransformation<Scalar>& expansion() const { return this->C; }


    /**
     *  @param sigma            alpha or beta
     * 
     *  @return the number of sigma electrons that these UHF model parameters describe, i.e. the number of occupied sigma-spin-orbitals
     */
    size_t numberOfElectrons(const Spin sigma) const {

        switch (sigma) {
        case Spin::alpha: {
            return this->N_alpha;
            break;
        }

        case Spin::beta: {
            return this->N_beta;
            break;
        }
        }
    }


    /**
     *  @return the total number of spin-orbitals that these UHF model parameters describe
     */
    size_t numberOfSpinOrbitals() const {

        return this->numberOfSpinOrbitals(Spin::alpha) + this->numberOfSpinOrbitals(Spin::beta);
    }


    /**
     *  @param sigma            alpha or beta
     * 
     *  @return the number of sigma spin-orbitals that these UHF model parameters describe
     */
    size_t numberOfSpinOrbitals(const Spin sigma) const {
        return this->expansion().component(sigma).numberOfOrbitals();
    }


    /**
     *  @param sigma            alpha or beta
     * 
     *  @return the orbital energies of the sigma-spin-orbitals
     */
    const VectorX<double>& orbitalEnergies(const Spin sigma) const {

        switch (sigma) {
        case Spin::alpha: {
            return this->orbital_energies_alpha;
            break;
        }

        case Spin::beta: {
            return this->orbital_energies_beta;
            break;
        }
        }
    }


    /**
     *  @return The implicit alpha and beta occupied-virtual orbital spaces that are associated to these UHF model parameters.
     */
    SpinResolvedOrbitalSpace orbitalSpace() const { return UHF<Scalar>::orbitalSpace(this->numberOfSpinOrbitals(Spin::alpha), this->numberOfElectrons(Spin::alpha),
                                                                                     this->numberOfSpinOrbitals(Spin::beta), this->numberOfElectrons(Spin::beta)); }


    /**
     *  @return all the spin-orbital energies, with the alpha spin-orbital energies appearing before the beta spin-orbital energies
     */
    VectorX<double> spinOrbitalEnergiesBlocked() const {

        GQCP::VectorX<double> total_orbital_energies {this->numberOfSpinOrbitals()};
        total_orbital_energies.head(this->numberOfSpinOrbitals(Spin::alpha)) = this->orbitalEnergies(Spin::alpha);
        total_orbital_energies.tail(this->numberOfSpinOrbitals(Spin::beta)) = this->orbitalEnergies(Spin::beta);

        return total_orbital_energies;
    }
};


}  // namespace QCModel
}  // namespace GQCP
