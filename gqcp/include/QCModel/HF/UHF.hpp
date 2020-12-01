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
#include "Operator/SecondQuantized/MixedUSQTwoElectronOperatorComponent.hpp"
#include "Operator/SecondQuantized/RSQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/USQOneElectronOperator.hpp"
#include "QCModel/HF/RHF.hpp"
#include "QCModel/HF/StabilityMatrices/UHFStabilityMatrices.hpp"
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
     * Construct a mixed-spin component of the spin-conserved stability matrix A'.
     * 
     *  @note The formula for this component is as follows:
     *      A'_I(sigma)A(sigma)J(sigma_bar)B(sigma_bar) = (A(sigma)I(sigma)|J(sigma_bar)B(sigma_bar)).
     * 
     *  @note Sigma_bar = alpha if sigma is beta and vice versa.
     * 
     *  @note The name `Spin-conserved`comes from the article "Constraints and stability in Hartree-Fock theory" by Seeger, R. and Pople J.A. (https://doi.org/10.1063/1.434318).
     * 
     *  @param usq_hamiltonian      The second quantized Hamiltonian, expressed in the orthonormal, 'unrestricted' spin orbital basis of the UHF MOs, which contains the necessary two-electron operators.
     *  @param sigma                The spin-sigma component. 
     *  @param sigma_bar            The spin component, opposite to sigma.
     * 
     *  @return The mixed spin sigma component of stability matrix A'.
     */
    GQCP::Matrix<Scalar> calculateMixedSpinConservedAComponent(const USQHamiltonian<Scalar>& usq_hamiltonian, const Spin sigma, const Spin sigma_bar) const {

        // Sigma and sigma_bar need to be different. If they are the same, the method throws an error.
        if (sigma == sigma_bar) {
            throw std::invalid_argument("QCModel::UHF<Scalar>.calculateMixedSpinConservedAComponent(const USQHamiltonian<Scalar>& usq_hamiltonian, const Spin sigma, const Spin sigma_bar): The spin 'sigma' and spin 'sigma_bar' arguments are not allowed to be the same.");
        }

        // Depending on whether we are making the aabb or bbaa A'-component, we need a different component of the orbital_space.
        const auto orbital_space = this->orbitalSpace();
        const auto& orbital_space_sigma = orbital_space.component(sigma);
        const auto& orbital_space_sigma_bar = orbital_space.component(sigma_bar);

        // Determine the number of occupied and virtual orbitals for both spin components.
        const auto& n_occ_sigma = orbital_space_sigma.numberOfOrbitals(OccupationType::k_occupied);
        const auto& n_virt_sigma = orbital_space_sigma.numberOfOrbitals(OccupationType::k_virtual);

        const auto& n_occ_sigma_bar = orbital_space_sigma_bar.numberOfOrbitals(OccupationType::k_occupied);
        const auto& n_virt_sigma_bar = orbital_space_sigma_bar.numberOfOrbitals(OccupationType::k_virtual);

        // The mixed alpha-beta two electron integrals are extracted from the Hamiltonian.
        const auto& g_aabb = usq_hamiltonian.twoElectron().alphaBeta().parameters();

        // The next step is to create the needed tensor slice.
        // Zero-initialize an occupied-virtual-occupied-virtual object of mixed spins.
        auto A_iajb_slice = orbital_space.template initializeMixedRepresentableObjectFor<Scalar>(OccupationType::k_occupied, sigma, OccupationType::k_virtual, sigma,
                                                                                                 OccupationType::k_occupied, sigma_bar, OccupationType::k_virtual, sigma_bar);
        for (const auto& i : orbital_space_sigma.indices(OccupationType::k_occupied)) {
            for (const auto& a : orbital_space_sigma.indices(OccupationType::k_virtual)) {
                for (const auto& j : orbital_space_sigma_bar.indices(OccupationType::k_occupied)) {
                    for (const auto& b : orbital_space_sigma_bar.indices(OccupationType::k_virtual)) {

                        // Depending on which spin component is alpha and which is beta, the indices need to be transposed correctly.
                        if (sigma == Spin::alpha && sigma_bar == Spin::beta) {
                            A_iajb_slice(i, a, j, b) = g_aabb(a, i, j, b);
                        } else {
                            A_iajb_slice(i, a, j, b) = g_aabb(j, b, a, i);
                        }
                    }
                }
            }
        }

        // Turn the ImplicitRankFourTensorSlice in an actual Tensor.
        auto A_iajb = A_iajb_slice.asTensor();

        // Finally, reshape the tensor to a matrix.
        const auto mixed_spin_conserved_A_component = A_iajb.reshape(n_occ_sigma * n_virt_sigma, n_occ_sigma_bar * n_virt_sigma_bar);

        return mixed_spin_conserved_A_component;
    }


    /**
     * Construct a mixed-spin component of the spin-conserved stability matrix B'.
     * 
     *  @note The formula for this component is as follows:
     *      B'_I(sigma)A(sigma)J(sigma_bar)B(sigma_bar) = (A(sigma)I(sigma)|B(sigma_bar)J(sigma_bar)).
     * 
     *  @note Sigma_bar = alpha if sigma is beta and vice versa.
     * 
     *  @note The name `Spin-conserved`comes from the article "Constraints and stability in Hartree-Fock theory" by Seeger, R. and Pople J.A. (https://doi.org/10.1063/1.434318).
     * 
     *  @param usq_hamiltonian      The second quantized Hamiltonian, expressed in the orthonormal, 'unrestricted' spin orbital basis of the UHF MOs, which contains the necessary two-electron operators.
     *  @param sigma                The spin-sigma component. 
     *  @param sigma_bar            The spin component, opposite to sigma.
     * 
     *  @return The mixed spin sigma component of stability matrix B'.
     */
    GQCP::Matrix<Scalar> calculateMixedSpinConservedBComponent(const USQHamiltonian<Scalar>& usq_hamiltonian, const Spin sigma, const Spin sigma_bar) const {

        // Sigma and sigma_bar need to be different. If they are the same, the method throws an error.
        if (sigma == sigma_bar) {
            throw std::invalid_argument("QCModel::UHF<Scalar>.calculateMixedSpinConservedBComponent(const USQHamiltonian<Scalar>& usq_hamiltonian, const Spin sigma, const Spin sigma_bar): The spin 'sigma' and spin 'sigma_bar' arguments are not allowed to be the same.");
        }

        // Depending on whether we are making the aabb or bbaa B'-component, we need a different component of the orbital_space.
        const auto orbital_space = this->orbitalSpace();
        const auto& orbital_space_sigma = orbital_space.component(sigma);
        const auto& orbital_space_sigma_bar = orbital_space.component(sigma_bar);

        // Determine the number of occupied and virtual orbitals for both spin components.
        const auto& n_occ_sigma = orbital_space_sigma.numberOfOrbitals(OccupationType::k_occupied);
        const auto& n_virt_sigma = orbital_space_sigma.numberOfOrbitals(OccupationType::k_virtual);

        const auto& n_occ_sigma_bar = orbital_space_sigma_bar.numberOfOrbitals(OccupationType::k_occupied);
        const auto& n_virt_sigma_bar = orbital_space_sigma_bar.numberOfOrbitals(OccupationType::k_virtual);

        // The mixed alpha-beta two electron integrals are extracted from the Hamiltonian.
        const auto& g_aabb = usq_hamiltonian.twoElectron().alphaBeta().parameters();

        // The next step is to create the needed tensor slice.
        // Zero-initialize an occupied-virtual-occupied-virtual object of mixed spins.
        auto B_iajb_slice = orbital_space.template initializeMixedRepresentableObjectFor<Scalar>(OccupationType::k_occupied, sigma, OccupationType::k_virtual, sigma,
                                                                                                 OccupationType::k_occupied, sigma_bar, OccupationType::k_virtual, sigma_bar);
        for (const auto& i : orbital_space_sigma.indices(OccupationType::k_occupied)) {
            for (const auto& a : orbital_space_sigma.indices(OccupationType::k_virtual)) {
                for (const auto& j : orbital_space_sigma_bar.indices(OccupationType::k_occupied)) {
                    for (const auto& b : orbital_space_sigma_bar.indices(OccupationType::k_virtual)) {

                        // Depending on which spin component is alpha and which is beta, the indices need to be transposed correctly.
                        if (sigma == Spin::alpha && sigma_bar == Spin::beta) {
                            B_iajb_slice(i, a, j, b) = g_aabb(a, i, b, j);
                        } else {
                            B_iajb_slice(i, a, j, b) = g_aabb(b, j, a, i);
                        }
                    }
                }
            }
        }

        // Turn the ImplicitRankFourTensorSlice in an actual Tensor.
        auto B_iajb = B_iajb_slice.asTensor();

        // Finally, reshape the tensor to a matrix.
        const auto mixed_spin_conserved_B_component = B_iajb.reshape(n_occ_sigma * n_virt_sigma, n_occ_sigma_bar * n_virt_sigma_bar);

        return mixed_spin_conserved_B_component;
    }


    /**
     * Construct a pure-spin component of the spin-conserved stability matrix A'.
     * 
     *  @note The formula for this component is as follows:
     *      A'_I(sigma)A(sigma)J(sigma)B(sigma) = \delta_I(sigma)J(sigma) * F_B(sigma)A(sigma) - \delta_B(sigma)A(sigma) * F_I(sigma)J(sigma) + (A(sigma)I(sigma)||J(sigma)B(sigma)).
     * 
     *  @note The name `Spin-conserved`comes from the article "Constraints and stability in Hartree-Fock theory" by Seeger, R. and Pople J.A. (https://doi.org/10.1063/1.434318).
     * 
     *  @param usq_hamiltonian      The second quantized Hamiltonian, expressed in the orthonormal, 'unrestricted' spin orbital basis of the UHF MOs, which contains the necessary two-electron operators.
     *  @param sigma                The spin-sigma component. This method creates a pure spin component, so all sigma's in the formula are either alpha or beta.
     * 
     *  @return The pure spin sigma component of stability matrix A'.
     */
    GQCP::Matrix<Scalar> calculatePureSpinConservedAComponent(const USQHamiltonian<Scalar>& usq_hamiltonian, const Spin sigma) const {

        // Depending on whether we are making the alpha or beta A'-component, we need a different component of the orbital_space.
        const auto orbital_space_sigma = this->orbitalSpace().component(sigma);

        // Determine the number of occupied and virtual orbitals.
        const auto& n_occ = orbital_space_sigma.numberOfOrbitals(OccupationType::k_occupied);
        const auto& n_virt = orbital_space_sigma.numberOfOrbitals(OccupationType::k_virtual);

        // The pure alpha-alpha two electron integrals are extracted from the Hamiltonian.
        // We need the anti-symmetrized tensor: (AI||JB) = (AI|JB) - (AB|JI). This is obtained by the `.antisymmetrized()` method.
        const auto g_aaaa = usq_hamiltonian.twoElectron().alphaAlpha().antisymmetrized().parameters();

        // The elements F_BA and F_IJ are the eigenvalues of the one-electron Fock operator.
        // The excitationEnergies API can be used to find these values.
        const auto F_values = this->excitationEnergies().component(sigma);

        // The next step is to create the needed tensor slice.
        // Zero-initialize an occupied-virtual-occupied-virtual object.
        auto A_iajb_slice = orbital_space_sigma.template initializeRepresentableObjectFor<Scalar>(OccupationType::k_occupied, OccupationType::k_virtual, OccupationType::k_occupied, OccupationType::k_virtual);
        for (const auto& i : orbital_space_sigma.indices(OccupationType::k_occupied)) {
            for (const auto& a : orbital_space_sigma.indices(OccupationType::k_virtual)) {
                for (const auto& j : orbital_space_sigma.indices(OccupationType::k_occupied)) {
                    for (const auto& b : orbital_space_sigma.indices(OccupationType::k_virtual)) {
                        A_iajb_slice(i, a, j, b) = g_aaaa(a, i, j, b);
                    }
                }
            }
        }

        // Turn the ImplicitRankFourTensorSlice in an actual Tensor.
        auto A_iajb = A_iajb_slice.asTensor();

        // Add the previously calculated F values on the correct positions.
        for (int a = 0; a < n_virt; a++) {
            for (int i = 0; i < n_occ; i++) {
                A_iajb(i, a, i, a) += F_values(a, i);
            }
        }

        // Finally, reshape the tensor to a matrix.
        const auto pure_spin_conserved_A_component = A_iajb.reshape(n_occ * n_virt, n_occ * n_virt);

        return pure_spin_conserved_A_component;
    }


    /**
     * Construct a pure-spin component of the spin-conserved stability matrix B'.
     * 
     *  @note The formula for this component is as follows:
     *      B'_I(sigma)A(sigma)J(sigma)B(sigma) = (A(sigma)I(sigma)||B(sigma)J(sigma)).
     * 
     *  @note The name `Spin-conserved`comes from the article "Constraints and stability in Hartree-Fock theory" by Seeger, R. and Pople J.A. (https://doi.org/10.1063/1.434318).
     * 
     *  @param usq_hamiltonian      The second quantized Hamiltonian, expressed in the orthonormal, 'unrestricted' spin orbital basis of the UHF MOs, which contains the necessary two-electron operators.
     *  @param sigma                The spin-sigma component. This method creates a pure spin component, so all sigma's in the formula are either alpha or beta.
     * 
     *  @return The pure spin sigma component of stability matrix B'.
     */
    GQCP::Matrix<Scalar> calculatePureSpinConservedBComponent(const USQHamiltonian<Scalar>& usq_hamiltonian, const Spin sigma) const {

        // Depending on whether we are making the alpha or beta B'-component, we need a different component of the orbital_space.
        const auto orbital_space_sigma = this->orbitalSpace().component(sigma);

        // Determine the number of occupied and virtual orbitals.
        const auto& n_occ = orbital_space_sigma.numberOfOrbitals(OccupationType::k_occupied);
        const auto& n_virt = orbital_space_sigma.numberOfOrbitals(OccupationType::k_virtual);

        // The pure beta-beta two electron integrals are extracted from the Hamiltonian.
        // We need the anti-symmetrized tensor: (AI||JB) = (AI|JB) - (AB|JI). This is obtained by the `.antisymmetrized()` method.
        const auto g_bbbb = usq_hamiltonian.twoElectron().betaBeta().antisymmetrized().parameters();

        // The next step is to create the needed tensor slice.
        // Zero-initialize an occupied-virtual-occupied-virtual object.
        auto B_iajb_slice = orbital_space_sigma.template initializeRepresentableObjectFor<Scalar>(OccupationType::k_occupied, OccupationType::k_virtual, OccupationType::k_occupied, OccupationType::k_virtual);
        for (const auto& i : orbital_space_sigma.indices(OccupationType::k_occupied)) {
            for (const auto& a : orbital_space_sigma.indices(OccupationType::k_virtual)) {
                for (const auto& j : orbital_space_sigma.indices(OccupationType::k_occupied)) {
                    for (const auto& b : orbital_space_sigma.indices(OccupationType::k_virtual)) {
                        B_iajb_slice(i, a, j, b) = g_bbbb(a, i, b, j);
                    }
                }
            }
        }

        // Turn the ImplicitRankFourTensorSlice in an actual Tensor.
        auto B_iajb = B_iajb_slice.asTensor();

        // Finally, reshape the tensor to a matrix.
        const auto pure_spin_conserved_B_component = B_iajb.reshape(n_occ * n_virt, n_occ * n_virt);

        return pure_spin_conserved_B_component;
    }


    /**
     * Construct the complete spin-conserved stability matrix A'.
     * 
     *  @note The formula for this component is as follows:
     *      A' = (A'_aaaa   A'_aabb)
     *           (A'_bbaa   A'_bbbb).
     * 
     *  @note The name `Spin-conserved`comes from the article "Constraints and stability in Hartree-Fock theory" by Seeger, R. and Pople J.A. (https://doi.org/10.1063/1.434318).
     * 
     *  @param usq_hamiltonian      The second quantized Hamiltonian, expressed in the orthonormal, 'unrestricted' spin orbital basis of the UHF MOs, which contains the necessary two-electron operators.
     * 
     *  @return The complete stability matrix A'.
     */
    GQCP::MatrixX<Scalar> calculateSpinConservedA(const USQHamiltonian<Scalar>& usq_hamiltonian) const {

        // Calculate the four different A' components.
        const auto A_aaaa = this->calculatePureSpinConservedAComponent(usq_hamiltonian, Spin::alpha);               // Dimension = (n_occ_a * n_virt_a, n_occ_a * n_virt_a).
        const auto A_bbbb = this->calculatePureSpinConservedAComponent(usq_hamiltonian, Spin::beta);                // Dimension = (n_occ_b * n_virt_b, n_occ_b * n_virt_b).
        const auto A_aabb = this->calculateMixedSpinConservedAComponent(usq_hamiltonian, Spin::alpha, Spin::beta);  // Dimension = (n_occ_a * n_virt_a, n_occ_b * n_virt_b).
        const auto A_bbaa = this->calculateMixedSpinConservedAComponent(usq_hamiltonian, Spin::beta, Spin::alpha);  // Dimension = (n_occ_b * n_virt_b, n_occ_a * n_virt_a).

        // Determine the total matrix dimension and initialize the total matrix.
        const auto orbital_space = this->orbitalSpace();
        const auto& n_occ_a = orbital_space.component(Spin::alpha).numberOfOrbitals(OccupationType::k_occupied);
        const auto& n_virt_a = orbital_space.component(Spin::alpha).numberOfOrbitals(OccupationType::k_virtual);
        const auto& n_occ_b = orbital_space.component(Spin::beta).numberOfOrbitals(OccupationType::k_occupied);
        const auto& n_virt_b = orbital_space.component(Spin::beta).numberOfOrbitals(OccupationType::k_virtual);

        const auto dimension = (n_occ_a * n_virt_a) + (n_occ_b * n_virt_b);
        GQCP::MatrixX<Scalar> total_A {dimension, dimension};

        // Place the components on the correct positions in the total matrix.
        total_A.topLeftCorner(n_occ_a * n_virt_a, n_occ_a * n_virt_a) = A_aaaa;
        total_A.topRightCorner(n_occ_a * n_virt_a, n_occ_b * n_virt_b) = A_aabb;
        total_A.bottomLeftCorner(n_occ_b * n_virt_b, n_occ_a * n_virt_a) = A_bbaa;
        total_A.bottomRightCorner(n_occ_b * n_virt_b, n_occ_b * n_virt_b) = A_bbbb;

        return total_A;
    }


    /**
     * Construct the complete spin-conserved stability matrix B'.
     * 
     *  @note The formula for this component is as follows:
     *      B' = (B'_aaaa   B'_aabb)
     *           (B'_bbaa   B'_bbbb).
     * 
     *  @note The name `Spin-conserved`comes from the article "Constraints and stability in Hartree-Fock theory" by Seeger, R. and Pople J.A. (https://doi.org/10.1063/1.434318).
     * 
     *  @param usq_hamiltonian      The second quantized Hamiltonian, expressed in the orthonormal, 'unrestricted' spin orbital basis of the UHF MOs, which contains the necessary two-electron operators.
     * 
     *  @return The complete stability matrix B'.
     */
    GQCP::MatrixX<Scalar> calculateSpinConservedB(const USQHamiltonian<Scalar>& usq_hamiltonian) const {

        // Calculate the four different A' components.
        const auto B_aaaa = this->calculatePureSpinConservedBComponent(usq_hamiltonian, Spin::alpha);               // Dimension = (n_occ_a * n_virt_a, n_occ_a * n_virt_a).
        const auto B_bbbb = this->calculatePureSpinConservedBComponent(usq_hamiltonian, Spin::beta);                // Dimension = (n_occ_b * n_virt_b, n_occ_b * n_virt_b).
        const auto B_aabb = this->calculateMixedSpinConservedBComponent(usq_hamiltonian, Spin::alpha, Spin::beta);  // Dimension = (n_occ_a * n_virt_a, n_occ_b * n_virt_b).
        const auto B_bbaa = this->calculateMixedSpinConservedBComponent(usq_hamiltonian, Spin::beta, Spin::alpha);  // Dimension = (n_occ_b * n_virt_b, n_occ_a * n_virt_a).

        // Determine the total matrix dimension and initialize the total matrix.
        const auto orbital_space = this->orbitalSpace();
        const auto& n_occ_a = orbital_space.component(Spin::alpha).numberOfOrbitals(OccupationType::k_occupied);
        const auto& n_virt_a = orbital_space.component(Spin::alpha).numberOfOrbitals(OccupationType::k_virtual);
        const auto& n_occ_b = orbital_space.component(Spin::beta).numberOfOrbitals(OccupationType::k_occupied);
        const auto& n_virt_b = orbital_space.component(Spin::beta).numberOfOrbitals(OccupationType::k_virtual);

        const auto dimension = (n_occ_a * n_virt_a) + (n_occ_b * n_virt_b);
        GQCP::MatrixX<Scalar> total_B {dimension, dimension};

        // Place the components on the correct positions in the total matrix.
        total_B.topLeftCorner(n_occ_a * n_virt_a, n_occ_a * n_virt_a) = B_aaaa;
        total_B.topRightCorner(n_occ_a * n_virt_a, n_occ_b * n_virt_b) = B_aabb;
        total_B.bottomLeftCorner(n_occ_b * n_virt_b, n_occ_a * n_virt_a) = B_bbaa;
        total_B.bottomRightCorner(n_occ_b * n_virt_b, n_occ_b * n_virt_b) = B_bbbb;

        return total_B;
    }


    /**
     * Construct a component of the spin-unconserved stability matrix A''.
     * 
     *  @note The formula for this component is as follows:
     *      A''_I(sigma)A(sigma_bar)J(sigma)B(sigma_bar) = \delta_I(sigma)J(sigma) * F_B(sigma_bar)A(sigma_bar) - \delta_B(sigma_bar)A(sigma_bar) * F_I(sigma)J(sigma) - (A(sigma_bar)B(sigma_bar)|J(sigma)I(sigma)).
     * 
     *  @note Sigma_bar = alpha if sigma is beta and vice versa.
     * 
     *  @note The name `Spin-unconserved`comes from the article "Constraints and stability in Hartree-Fock theory" by Seeger, R. and Pople J.A. (https://doi.org/10.1063/1.434318).
     * 
     *  @param usq_hamiltonian      The second quantized Hamiltonian, expressed in the orthonormal, 'unrestricted' spin orbital basis of the UHF MOs, which contains the necessary two-electron operators.
     *  @param sigma                The spin-sigma component. 
     *  @param sigma_bar            The spin component, opposite to sigma.
     * 
     *  @return A component of stability matrix A''.
     */
    GQCP::Matrix<Scalar> calculateSpinUnconservedAComponent(const USQHamiltonian<Scalar>& usq_hamiltonian, const Spin sigma, const Spin sigma_bar) const {

        // Sigma and sigma_bar need to be different. If they are the same, the method throws an error.
        if (sigma == sigma_bar) {
            throw std::invalid_argument("QCModel::UHF<Scalar>.calculateMixedSpinConservedAComponent(const USQHamiltonian<Scalar>& usq_hamiltonian, const Spin sigma, const Spin sigma_bar): The spin 'sigma' and spin 'sigma_bar' arguments are not allowed to be the same.");
        }

        // Depending on whether we're making the abab or baba A''-component, we need a different component of the orbital_space.
        const auto orbital_space = this->orbitalSpace();
        const auto& orbital_space_sigma = orbital_space.component(sigma);
        const auto& orbital_space_sigma_bar = orbital_space.component(sigma_bar);

        // Determine the number of occupied and virtual orbitals.
        const auto& n_occ_sigma = orbital_space_sigma.numberOfOrbitals(OccupationType::k_occupied);
        const auto& n_virt_sigma = orbital_space_sigma.numberOfOrbitals(OccupationType::k_virtual);

        const auto& n_occ_sigma_bar = orbital_space_sigma_bar.numberOfOrbitals(OccupationType::k_occupied);
        const auto& n_virt_sigma_bar = orbital_space_sigma_bar.numberOfOrbitals(OccupationType::k_virtual);

        /// The mixed alpha-beta two electron integrals are extracted from the Hamiltonian.
        const auto& g_aabb = usq_hamiltonian.twoElectron().alphaBeta().parameters();

        // The next step is to create the needed tensor slice.
        // Zero-initialize an occupied-virtual-occupied-virtual object of mixed spins.
        auto A_iajb_slice = orbital_space.template initializeMixedRepresentableObjectFor<Scalar>(OccupationType::k_occupied, sigma_bar, OccupationType::k_virtual, sigma,
                                                                                                 OccupationType::k_occupied, sigma_bar, OccupationType::k_virtual, sigma);
        for (const auto& i : orbital_space_sigma_bar.indices(OccupationType::k_occupied)) {
            for (const auto& a : orbital_space_sigma.indices(OccupationType::k_virtual)) {
                for (const auto& j : orbital_space_sigma_bar.indices(OccupationType::k_occupied)) {
                    for (const auto& b : orbital_space_sigma.indices(OccupationType::k_virtual)) {

                        // Depending on which spin component is alpha and which is beta, the indices need to be transposed correctly.
                        if (sigma == Spin::alpha && sigma_bar == Spin::beta) {
                            A_iajb_slice(i, a, j, b) = -g_aabb(a, b, j, i);
                        } else {
                            A_iajb_slice(i, a, j, b) = -g_aabb(j, i, a, b);
                        }
                    }
                }
            }
        }
        // Turn the ImplicitRankFourTensorSlice in an actual Tensor.
        auto A_iajb = A_iajb_slice.asTensor();

        // The elements F_BA and F_IJ are the eigenvalues of the one-electron Fock operator.
        // We have to construct a mixed matrix in this particular case.
        GQCP::MatrixX<Scalar> F_values {n_virt_sigma, n_occ_sigma_bar};
        const auto virtual_energies = this->virtualOrbitalEnergies().component(sigma);
        const auto occupied_energies = this->occupiedOrbitalEnergies().component(sigma_bar);

        for (int a = 0; a < n_virt_sigma; a++) {
            for (int i = 0; i < n_occ_sigma_bar; i++) {
                F_values(a, i) = virtual_energies[a] - occupied_energies[i];
            }
        }

        // Add the previously calculated F values on the correct positions.
        for (int a = 0; a < n_virt_sigma; a++) {
            for (int i = 0; i < n_occ_sigma_bar; i++) {
                A_iajb(i, a, i, a) += F_values(a, i);
            }
        }

        // Finally, reshape the tensor to a matrix.
        const auto spin_unconserved_A_component = A_iajb.reshape(n_occ_sigma_bar * n_virt_sigma, n_occ_sigma_bar * n_virt_sigma);

        return spin_unconserved_A_component;
    }


    /**
     * Construct a component of the spin-unconserved stability matrix B''.
     * 
     *  @note The formula for this component is as follows:
     *      B''_I(sigma)A(sigma_bar)J(sigma_bar)B(sigma) = - (A(sigma_bar)J(sigma_bar)|B(sigma)I(sigma)).
     * 
     *  @note Sigma_bar = alpha if sigma is beta and vice versa.
     * 
     *  @note The name `Spin-unconserved`comes from the article "Constraints and stability in Hartree-Fock theory" by Seeger, R. and Pople J.A. (https://doi.org/10.1063/1.434318).
     * 
     *  @param usq_hamiltonian      The second quantized Hamiltonian, expressed in the orthonormal, 'unrestricted' spin orbital basis of the UHF MOs, which contains the necessary two-electron operators.
     *  @param sigma                The spin-sigma component. 
     *  @param sigma_bar            The spin component, opposite to sigma.
     * 
     *  @return A component of stability matrix B''.
     */
    GQCP::Matrix<Scalar> calculateSpinUnconservedBComponent(const USQHamiltonian<Scalar>& usq_hamiltonian, const Spin sigma, const Spin sigma_bar) const {

        // Sigma and sigma_bar need to be different. If they are the same, the method throws an error.
        if (sigma == sigma_bar) {
            throw std::invalid_argument("QCModel::UHF<Scalar>.calculateMixedSpinConservedAComponent(const USQHamiltonian<Scalar>& usq_hamiltonian, const Spin sigma, const Spin sigma_bar): The spin 'sigma' and spin 'sigma_bar' arguments are not allowed to be the same.");
        }

        // Depending on whether we're making the abba or baab B''-component, we need a different component of the orbital_space.
        const auto orbital_space = this->orbitalSpace();
        const auto& orbital_space_sigma = orbital_space.component(sigma);
        const auto& orbital_space_sigma_bar = orbital_space.component(sigma_bar);

        // Determine the number of occupied and virtual orbitals.
        const auto& n_occ_sigma = orbital_space_sigma.numberOfOrbitals(OccupationType::k_occupied);
        const auto& n_virt_sigma = orbital_space_sigma.numberOfOrbitals(OccupationType::k_virtual);

        const auto& n_occ_sigma_bar = orbital_space_sigma_bar.numberOfOrbitals(OccupationType::k_occupied);
        const auto& n_virt_sigma_bar = orbital_space_sigma_bar.numberOfOrbitals(OccupationType::k_virtual);

        // The mixed alpha-beta two electron integrals are extracted from the Hamiltonian.
        const auto& g_aabb = usq_hamiltonian.twoElectron().alphaBeta().parameters();

        // The next step is to create the needed tensor slice.
        // Zero-initialize an occupied-virtual-occupied-virtual object of mixed spins.
        auto B_iajb_slice = orbital_space.template initializeMixedRepresentableObjectFor<Scalar>(OccupationType::k_occupied, sigma, OccupationType::k_virtual, sigma_bar,
                                                                                                 OccupationType::k_occupied, sigma_bar, OccupationType::k_virtual, sigma);
        for (const auto& i : orbital_space_sigma.indices(OccupationType::k_occupied)) {
            for (const auto& a : orbital_space_sigma_bar.indices(OccupationType::k_virtual)) {
                for (const auto& j : orbital_space_sigma_bar.indices(OccupationType::k_occupied)) {
                    for (const auto& b : orbital_space_sigma.indices(OccupationType::k_virtual)) {

                        // Depending on which spin component is alpha and which is beta, the indices need to be transposed correctly.
                        if (sigma == Spin::alpha && sigma_bar == Spin::beta) {
                            B_iajb_slice(i, a, j, b) = -g_aabb(b, i, a, j);
                        } else {
                            B_iajb_slice(i, a, j, b) = -g_aabb(a, j, b, i);
                        }
                    }
                }
            }
        }

        // Turn the ImplicitRankFourTensorSlice in an actual Tensor
        auto B_iajb = B_iajb_slice.asTensor();

        // Finally, reshape the tensor to a matrix.

        const auto spin_unconserved_B_component = B_iajb.reshape(n_occ_sigma_bar * n_virt_sigma, n_occ_sigma * n_virt_sigma_bar);

        return spin_unconserved_B_component;
    }


    /**
     * Construct the complete spin-unconserved stability matrix A''.
     * 
     *  @note The formula for this component is as follows:
     *      A'' = (A''_abab      0    )
     *            (   0       A''_baba).
     * 
     *  @note The name `Spin-unconserved`comes from the article "Constraints and stability in Hartree-Fock theory" by Seeger, R. and Pople J.A. (https://doi.org/10.1063/1.434318).
     * 
     *  @param usq_hamiltonian      The second quantized Hamiltonian, expressed in the orthonormal, 'unrestricted' spin orbital basis of the UHF MOs, which contains the necessary two-electron operators.
     * 
     *  @return The complete stability matrix A''.
     */
    GQCP::MatrixX<Scalar> calculateSpinUnconservedA(const USQHamiltonian<Scalar>& usq_hamiltonian) const {

        // Calculate the two different A'' components.
        const auto A_abab = this->calculateSpinUnconservedAComponent(usq_hamiltonian, Spin::alpha, Spin::beta);  // Dimension = (n_occ_b * n_virt_a, n_occ_b * n_virt_a).
        const auto A_baba = this->calculateSpinUnconservedAComponent(usq_hamiltonian, Spin::beta, Spin::alpha);  // Dimension = (n_occ_a * n_virt_b, n_occ_a * n_virt_b).

        // Determine the total matrix dimension and initialize the total matrix.
        const auto orbital_space = this->orbitalSpace();
        const auto& n_occ_a = orbital_space.component(Spin::alpha).numberOfOrbitals(OccupationType::k_occupied);
        const auto& n_virt_a = orbital_space.component(Spin::alpha).numberOfOrbitals(OccupationType::k_virtual);
        const auto& n_occ_b = orbital_space.component(Spin::beta).numberOfOrbitals(OccupationType::k_occupied);
        const auto& n_virt_b = orbital_space.component(Spin::beta).numberOfOrbitals(OccupationType::k_virtual);

        const auto dimension = (n_occ_b * n_virt_a) + (n_occ_a * n_virt_b);
        GQCP::MatrixX<Scalar> total_A {dimension, dimension};

        // Initialize the zero blocks of the total matrix.
        const auto zero_1 = GQCP::MatrixX<Scalar>::Zero(n_virt_a * n_occ_b, n_occ_a * n_virt_b);
        const auto zero_2 = GQCP::MatrixX<Scalar>::Zero(n_occ_a * n_virt_b, n_virt_a * n_occ_b);

        // Place the components on the correct positions in the total matrix.
        total_A.topLeftCorner(n_occ_b * n_virt_a, n_occ_b * n_virt_a) = A_abab;
        total_A.topRightCorner(n_virt_a * n_occ_b, n_occ_a * n_virt_b) = zero_1;
        total_A.bottomLeftCorner(n_occ_a * n_virt_b, n_occ_b * n_virt_a) = zero_2;
        total_A.bottomRightCorner(n_occ_a * n_virt_b, n_occ_a * n_virt_b) = A_baba;

        return total_A;
    }


    /**
     * Construct the complete spin-unconserved stability matrix B''.
     * 
     *  @note The formula for this component is as follows:
     *      B'' = (   0     B''_abba)
     *            (B''_baab    0    ).
     * 
     *  @note The name `Spin-unconserved`comes from the article "Constraints and stability in Hartree-Fock theory" by Seeger, R. and Pople J.A. (https://doi.org/10.1063/1.434318).
     * 
     *  @param usq_hamiltonian      The second quantized Hamiltonian, expressed in the orthonormal, 'unrestricted' spin orbital basis of the UHF MOs, which contains the necessary two-electron operators.
     * 
     *  @return The complete stability matrix B''.
     */
    GQCP::MatrixX<Scalar> calculateSpinUnconservedB(const USQHamiltonian<Scalar>& usq_hamiltonian) const {

        // Calculate the two different B'' components.
        const auto B_abba = this->calculateSpinUnconservedBComponent(usq_hamiltonian, Spin::alpha, Spin::beta);  // Dimension = (n_occ_b * n_virt_a, n_occ_a * n_virt_b).
        const auto B_baab = this->calculateSpinUnconservedBComponent(usq_hamiltonian, Spin::beta, Spin::alpha);  // Dimension = (n_occ_a * n_virt_b, n_occ_b * n_virt_a).

        // Determine the total matrix dimension and initialize the total matrix.
        const auto orbital_space = this->orbitalSpace();
        const auto& n_occ_a = orbital_space.component(Spin::alpha).numberOfOrbitals(OccupationType::k_occupied);
        const auto& n_virt_a = orbital_space.component(Spin::alpha).numberOfOrbitals(OccupationType::k_virtual);
        const auto& n_occ_b = orbital_space.component(Spin::beta).numberOfOrbitals(OccupationType::k_occupied);
        const auto& n_virt_b = orbital_space.component(Spin::beta).numberOfOrbitals(OccupationType::k_virtual);

        const auto dimension = (n_occ_b * n_virt_a) + (n_occ_a * n_virt_b);
        GQCP::MatrixX<Scalar> total_B {dimension, dimension};

        // Initialize the zero blocks of the total matrix.
        const auto zero_1 = GQCP::MatrixX<Scalar>::Zero(n_virt_a * n_occ_b, n_occ_b * n_virt_a);
        const auto zero_2 = GQCP::MatrixX<Scalar>::Zero(n_occ_a * n_virt_b, n_virt_b * n_occ_a);

        // Place the components on the correct positions in the total matrix.
        total_B.topLeftCorner(n_occ_b * n_virt_a, n_occ_b * n_virt_a) = zero_1;
        total_B.topRightCorner(n_virt_a * n_occ_b, n_occ_a * n_virt_b) = B_abba;
        total_B.bottomLeftCorner(n_occ_a * n_virt_b, n_occ_b * n_virt_a) = B_baab;
        total_B.bottomRightCorner(n_occ_a * n_virt_b, n_occ_a * n_virt_b) = zero_2;

        return total_B;
    }


    /**
     *  Calculate the UHF stability matrices and return them.
     * 
     *  @param usq_hamiltonian      The second quantized Hamiltonian, expressed in the orthonormal, 'unrestricted' spin orbital basis of the UHF MOs, which contains the necessary two-electron operators.
     *
     *  @return The UHF stability matrices.
     */
    UHFStabilityMatrices<Scalar> calculateStabilityMatrices(const USQHamiltonian<Scalar>& usq_hamiltonian) const {
        return UHFStabilityMatrices<Scalar> {this->calculateSpinConservedA(usq_hamiltonian), this->calculateSpinConservedB(usq_hamiltonian),
                                             this->calculateSpinUnconservedA(usq_hamiltonian), this->calculateSpinUnconservedB(usq_hamiltonian)};
    }


    /**
     *  @param sigma                The spin-sigma component. 
     * 
     *  @return A matrix containing all the possible excitation energies of the wavefunction model, belonging to a certain spin component. 
     * 
     *  @note       The rows are determined by the number of virtual sigma orbitals, the columns by the number of occupied sigma orbitals.
     * 
     */
    SpinResolved<GQCP::MatrixX<Scalar>> excitationEnergies() const {

        // Create the orbital space to determine the loops.
        const auto orbital_space = this->orbitalSpace();

        // Determine the number of occupied and virtual orbitals.
        const auto& n_occ_a = orbital_space.alpha().numberOfOrbitals(OccupationType::k_occupied);
        const auto& n_virt_a = orbital_space.alpha().numberOfOrbitals(OccupationType::k_virtual);
        const auto& n_occ_b = orbital_space.beta().numberOfOrbitals(OccupationType::k_occupied);
        const auto& n_virt_b = orbital_space.beta().numberOfOrbitals(OccupationType::k_virtual);

        // Calculate the occupied and virtual orbital energies.
        const auto occupied_energies = this->occupiedOrbitalEnergies();
        const auto virtual_energies = this->virtualOrbitalEnergies();

        // Create the F matrix.
        GQCP::MatrixX<Scalar> F_values_a {n_virt_a, n_occ_a};
        GQCP::MatrixX<Scalar> F_values_b {n_virt_b, n_occ_b};

        for (int a = 0; a < n_virt_a; a++) {
            for (int i = 0; i < n_occ_a; i++) {
                F_values_a(a, i) = virtual_energies.alpha()[a] - occupied_energies.alpha()[i];
            }
        }
        for (int a = 0; a < n_virt_b; a++) {
            for (int i = 0; i < n_occ_b; i++) {
                F_values_b(a, i) = virtual_energies.beta()[a] - occupied_energies.beta()[i];
            }
        }
        return SpinResolved<GQCP::MatrixX<Scalar>> {F_values_a, F_values_b};
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
     *  @param sigma            The alpha or beta spin component.
     * 
     *  @return The orbital energies belonging to the occupied orbitals of spin component sigma.
     */
    SpinResolved<std::vector<double>> occupiedOrbitalEnergies() const {

        // Determine the number of occupied orbitals of the spin components.
        const auto orbital_space = this->orbitalSpace();
        const auto& n_occ_a = orbital_space.alpha().numberOfOrbitals(OccupationType::k_occupied);
        const auto& n_occ_b = orbital_space.beta().numberOfOrbitals(OccupationType::k_occupied);


        std::vector<double> mo_energies_a;  // We use a std::vector in order to be able to slice the vector later on.
        std::vector<double> mo_energies_b;  // We use a std::vector in order to be able to slice the vector later on.

        for (int i = 0; i < this->numberOfSpinOrbitals(Spin::alpha); i++) {
            mo_energies_a.push_back(this->orbitalEnergies(Spin::alpha)[i]);
        }
        for (int i = 0; i < this->numberOfSpinOrbitals(Spin::beta); i++) {
            mo_energies_b.push_back(this->orbitalEnergies(Spin::beta)[i]);
        }

        // Add the values with indices greater than the occupied orbital indices, i.e. the virtual orbital indices, to the new vector.
        std::vector<double> mo_energies_occupied_a;
        std::vector<double> mo_energies_occupied_b;

        // Add the values with indices smaller than the occupied orbital indices, to the new vector.
        std::copy(mo_energies_a.begin(), mo_energies_a.begin() + n_occ_a, std::back_inserter(mo_energies_occupied_a));
        std::copy(mo_energies_b.begin(), mo_energies_b.begin() + n_occ_b, std::back_inserter(mo_energies_occupied_b));

        return SpinResolved<std::vector<double>> {mo_energies_occupied_a, mo_energies_occupied_b};
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
    SpinResolvedOrbitalSpace orbitalSpace() const { return UHF<Scalar>::orbitalSpace(this->numberOfSpinOrbitals(Spin::alpha), this->numberOfSpinOrbitals(Spin::beta),
                                                                                     this->numberOfElectrons(Spin::alpha), this->numberOfElectrons(Spin::beta)); }


    /**
     *  @return all the spin-orbital energies, with the alpha spin-orbital energies appearing before the beta spin-orbital energies
     */
    VectorX<double> spinOrbitalEnergiesBlocked() const {

        GQCP::VectorX<double> total_orbital_energies {this->numberOfSpinOrbitals()};
        total_orbital_energies.head(this->numberOfSpinOrbitals(Spin::alpha)) = this->orbitalEnergies(Spin::alpha);
        total_orbital_energies.tail(this->numberOfSpinOrbitals(Spin::beta)) = this->orbitalEnergies(Spin::beta);

        return total_orbital_energies;
    }


    /**
     *  @param sigma            The alpha or beta spin component.
     * 
     *  @return The orbital energies belonging to the virtual orbitals of spin component sigma.
     */
    SpinResolved<std::vector<double>> virtualOrbitalEnergies() const {

        // Determine the number of occupied orbitals of the spin components.
        const auto orbital_space = this->orbitalSpace();
        const auto& n_occ_a = orbital_space.alpha().numberOfOrbitals(OccupationType::k_occupied);
        const auto& n_occ_b = orbital_space.beta().numberOfOrbitals(OccupationType::k_occupied);


        std::vector<double> mo_energies_a;  // We use a std::vector in order to be able to slice the vector later on.
        std::vector<double> mo_energies_b;  // We use a std::vector in order to be able to slice the vector later on.

        for (int i = 0; i < this->numberOfSpinOrbitals(Spin::alpha); i++) {
            mo_energies_a.push_back(this->orbitalEnergies(Spin::alpha)[i]);
        }
        for (int i = 0; i < this->numberOfSpinOrbitals(Spin::beta); i++) {
            mo_energies_b.push_back(this->orbitalEnergies(Spin::beta)[i]);
        }

        // Add the values with indices greater than the occupied orbital indices, i.e. the virtual orbital indices, to the new vector.
        std::vector<double> mo_energies_virtual_a;
        std::vector<double> mo_energies_virtual_b;

        std::copy(mo_energies_a.begin() + n_occ_a, mo_energies_a.end(), std::back_inserter(mo_energies_virtual_a));
        std::copy(mo_energies_b.begin() + n_occ_b, mo_energies_b.end(), std::back_inserter(mo_energies_virtual_b));

        return SpinResolved<std::vector<double>> {mo_energies_virtual_a, mo_energies_virtual_b};
    }
};


}  // namespace QCModel
}  // namespace GQCP
