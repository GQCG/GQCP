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


#include "Basis/Transformations/GTransformationMatrix.hpp"
#include "DensityMatrix/G1DM.hpp"
#include "Operator/FirstQuantized/ElectronicSpinOperator.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCModel/HF/StabilityMatrices/GHFStabilityMatrices.hpp"
#include "Utilities/aliases.hpp"


namespace GQCP {
namespace QCModel {


/**
 *  The generalized Hartree-Fock wave function model.
 * 
 *  @tparam _Scalar             the type of scalar that is used for the expansion of the spinors in their underlying scalar basis
 */
template <typename _Scalar>
class GHF {
public:
    using Scalar = _Scalar;

private:
    size_t N;  // the number of electrons

    VectorX<double> orbital_energies;  // sorted in ascending energies
    GTransformationMatrix<Scalar> C;   // the coefficient matrix that expresses every spinor (as a column) in the underlying scalar bases

public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  The standard member-wise constructor
     * 
     *  @param N                    the number of electrons
     *  @param C                    the coefficient matrix that expresses every spinor (as a column) in the underlying scalar bases
     *  @param orbital_energies     the GHF MO energies
     */
    GHF(const size_t N, const VectorX<double>& orbital_energies, const GTransformationMatrix<Scalar>& C) :
        N {N},
        orbital_energies {orbital_energies},
        C {C} {}


    /**
     *  Default constructor setting everything to zero
     */
    GHF() :
        GHF(0, VectorX<double>::Zero(0), GTransformationMatrix<Scalar>::Zero(0, 0)) {}


    /*
     *  STATIC PUBLIC METHODS
     */

    /**
     *  @param F                the Fock matrix expressed in the scalar bases
     *  @param D                the GHF density matrix in the same scalar bases
     *  @param S                the (spin-blocked) overlap matrix of the scalar bases
     * 
     *  @return the GHF error matrix
     */
    static SquareMatrix<Scalar> calculateError(const SquareMatrix<Scalar>& F, const G1DM<Scalar>& P, const SquareMatrix<Scalar>& S) {
        return F * P * S - S * P * F;
    }


    /**
     *  @param P                the (spin-blocked) GHF density matrix in the scalar bases
     *  @param H_core           the (spin-blocked) core Hamiltonian expressed in the same scalar bases
     *  @param F                the (spin-blocked) Fock matrix in the same scalar bases
     *
     *  @return the GHF electronic energy
     */
    static double calculateElectronicEnergy(const G1DM<Scalar>& P, const ScalarGSQOneElectronOperator<Scalar>& H_core, const ScalarGSQOneElectronOperator<Scalar>& F) {

        // First, calculate the sum of H_core and F (this saves a contraction).
        const auto Z = H_core + F;

        // Convert the matrix Z to an GQCP::Tensor<double, 2> Z_tensor.
        // Einsum is only implemented for a tensor + a matrix, not for 2 matrices.
        Eigen::TensorMap<Eigen::Tensor<const Scalar, 2>> Z_t {Z.parameters().data(), P.rows(), P.cols()};
        Tensor<Scalar, 2> Z_tensor = Tensor<Scalar, 2>(Z_t);

        // To calculate the electronic energy, we must perform a double contraction (with prefactor 0.5).
        //      0.5 D(mu nu) Z(mu nu)
        // See knowdes: https://gqcg-res.github.io/knowdes/general-hartree-fock-theory.html

        Tensor<Scalar, 0> contraction = 0.5 * Z_tensor.template einsum<2>("ij, ij ->", P);

        // As the double contraction of two matrices is a scalar (a tensor of rank 0), we should access the value as (0).
        return contraction(0);
    }


    /**
     *  @param C            the coefficient matrix that expresses every spinor orbital (as a column) in the underlying scalar bases
     *  @param N            the number of electrons
     *
     *  @return the GHF 1-DM expressed in the underlying scalar basis
     */
    static G1DM<Scalar> calculateScalarBasis1DM(const GTransformationMatrix<double>& C, const size_t N) {

        const size_t M = C.dimension();
        const auto P_orthonormal = GHF<Scalar>::calculateOrthonormalBasis1DM(M, N);

        // Transform the 1-DM in an orthonormal basis to the underlying scalar basis
        return C.conjugate() * P_orthonormal * C.transpose();
    }


    /**
     *  Calculate the GHF direct (Coulomb) matrix.
     * 
     *  @param P                    the GHF density matrix expressed in the underlying scalar orbital basis
     *  @param sq_hamiltonian       the Hamiltonian expressed in the scalar (AO) basis, resulting from a quantization using a GSpinorBasis
     * 
     *  @return the GHF direct (Coulomb) matrix
     * 
     *  @note The scalar bases for the alpha- and beta-components must be the same.
     */
    static ScalarGSQOneElectronOperator<Scalar> calculateScalarBasisDirectMatrix(const G1DM<Scalar>& P, const GSQHamiltonian<Scalar>& sq_hamiltonian) {

        // To perform the contraction, we will first have to convert the density matrix into a Eigen::Tensor (since contractions are only implemented for tensors).
        // Since the two-electron integrals are spin-blocked (due to the nature of quantizing in a GSpinorBasis), the contractions must happen with a density matrix of the same dimension (M: the number of spinors). Therefore, we will construct a zero density matrix in which we only fill in one of the spin-blocks.
        const auto M = P.dimension();  // the total number of basis functions
        G1DM<Scalar> P_aa = G1DM<Scalar>::Zero(M);
        P_aa.topLeftCorner(M / 2, M / 2) = P.topLeftCorner(M / 2, M / 2);

        G1DM<Scalar> P_bb = G1DM<Scalar>::Zero(M);
        P_bb.bottomRightCorner(M / 2, M / 2) = P.bottomRightCorner(M / 2, M / 2);

        // Specify the contraction pairs for the direct contractions:
        //      P(rho lambda) (mu nu|rho lambda)
        // See knowdes: https://gqcg-res.github.io/knowdes/derivation-of-the-ghf-scf-equations-through-lagrange-multipliers.html
        const auto& g = sq_hamiltonian.twoElectron().parameters();
        const auto Ja = g.template einsum<2>("ijkl,kl->ij", P_aa).asMatrix();
        const auto Jb = g.template einsum<2>("ijkl,kl->ij", P_bb).asMatrix();


        return ScalarGSQOneElectronOperator<Scalar>(Ja + Jb);
    }


    /**
     *  Calculate the GHF exchange matrix.
     * 
     *  @param P                     the GHF density matrix expressed in the underlying scalar orbital bases
     *  @param sq_hamiltonian        the Hamiltonian expressed in the scalar (AO) basis, resulting from a quantization using a GSpinorBasis
     * 
     *  @return the UHF direct (Coulomb) matrix for spin sigma
     */
    static ScalarGSQOneElectronOperator<Scalar> calculateScalarBasisExchangeMatrix(const G1DM<Scalar>& P, const GSQHamiltonian<Scalar>& sq_hamiltonian) {

        // To perform the contraction, we will first have to convert the density matrix into a Eigen::Tensor (since contractions are only implemented for tensors).
        // Since the two-electron integrals are spin-blocked (due to the nature of quantizing in a GSpinorBasis), the contractions must happen with a density matrix of the same dimension (M: the number of spinors). Therefore, we will construct a zero density matrix in which we only fill in one of the spin-blocks.
        const auto M = P.dimension();  // the total number of basis functions

        G1DM<Scalar> P_aa = G1DM<Scalar>::Zero(M);
        P_aa.topLeftCorner(M / 2, M / 2) = P.topLeftCorner(M / 2, M / 2);

        G1DM<Scalar> P_ab = G1DM<Scalar>::Zero(M);
        P_ab.topRightCorner(M / 2, M / 2) = P.topRightCorner(M / 2, M / 2);

        G1DM<Scalar> P_ba = G1DM<Scalar>::Zero(M);
        P_ba.bottomLeftCorner(M / 2, M / 2) = P.bottomLeftCorner(M / 2, M / 2);

        G1DM<Scalar> P_bb = G1DM<Scalar>::Zero(M);
        P_bb.bottomRightCorner(M / 2, M / 2) = P.bottomRightCorner(M / 2, M / 2);


        // Specify the contraction pairs for the exchange contractions:
        //      P(lambda rho) (mu rho|lambda nu)
        const auto& g = sq_hamiltonian.twoElectron().parameters();
        const auto K_aa = g.template einsum<2>("ijkl,kj->il", P_aa).asMatrix();
        const auto K_ab = g.template einsum<2>("ijkl,kj->il", P_ba).asMatrix();
        const auto K_ba = g.template einsum<2>("ijkl,kj->il", P_ab).asMatrix();
        const auto K_bb = g.template einsum<2>("ijkl,kj->il", P_bb).asMatrix();

        // Each of the spin-blocks are calculated separately (while the other blocks are zero), so the total exchange matrix can be calculated as the sum of each part.
        return ScalarGSQOneElectronOperator<Scalar>(K_aa + K_ab + K_ba + K_bb);
    }


    /**
     *  Calculate the GHF Fock matrix F = H_core + G, in which G is a contraction of the density matrix and the two-electron integrals
     *
     *  @param P                    the (spin-blocked) GHF density matrix in the scalar bases
     *  @param sq_hamiltonian       the Hamiltonian expressed in the scalar (AO) basis, resulting from a quantization using a GSpinorBasis
     *
     *  @return the GHF Fock matrix expressed in the scalar basis
     */
    static ScalarGSQOneElectronOperator<Scalar> calculateScalarBasisFockMatrix(const G1DM<Scalar>& P, const GSQHamiltonian<Scalar>& sq_hamiltonian) {

        const auto& H_core = sq_hamiltonian.core();
        const auto J = QCModel::GHF<Scalar>::calculateScalarBasisDirectMatrix(P, sq_hamiltonian);
        const auto K = QCModel::GHF<Scalar>::calculateScalarBasisExchangeMatrix(P, sq_hamiltonian);

        return H_core + J - K;
    }


    /**
     *  @param M            the number of spinors
     *  @param N            the number of electrons
     *
     *  @return the GHF 1-DM expressed in an orthonormal spinor basis
     */
    static G1DM<Scalar>
    calculateOrthonormalBasis1DM(const size_t M, const size_t N) {

        // The 1-DM for GHF looks like (for M=5, N=3)
        //    1  0  0  0  0
        //    0  1  0  0  0
        //    0  0  1  0  0
        //    0  0  0  0  0
        //    0  0  0  0  0

        G1DM<Scalar> D_MO = G1DM<Scalar>::Zero(M);
        D_MO.topLeftCorner(N, N) = SquareMatrix<Scalar>::Identity(N);

        return D_MO;
    }


    /**
     *  @param M            the number of spinors
     *  @param N            the number of electrons
     * 
     *  @return the implicit (i.e. with ascending and contiguous orbital indices) occupied-virtual orbital space that corresponds to these GHF model parameters
     */
    static OrbitalSpace orbitalSpace(const size_t M, const size_t N) {

        return OrbitalSpace::Implicit({{OccupationType::k_occupied, N}, {OccupationType::k_virtual, M - N}});
    }


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @param spin_op                      the electronic spin operator
     *  @param S                            the (spin-blocked) overlap matrix of the underlying AO bases
     * 
     *  @return the expectation value of the electronic spin operator
     */
    Vector<complex, 3> calculateExpectationValueOf(const ElectronicSpinOperator& spin_op, const SquareMatrix<Scalar>& S) const {

        // Prepare some variables.
        const auto M = this->numberOfSpinors();
        const MatrixX<complex> C_alpha = this->coefficientMatrix().topRows(M / 2);
        const MatrixX<complex> C_beta = this->coefficientMatrix().bottomRows(M / 2);
        const SquareMatrix<complex> S_AO = S.topLeftCorner(M / 2, M / 2);  // assume equal for alpha and beta

        // Calculate overlaps between the alpha- and beta-spinors.
        const MatrixX<complex> overlap_aa = C_alpha.adjoint() * S_AO * C_alpha;
        const MatrixX<complex> overlap_ab = C_alpha.adjoint() * S_AO * C_beta;
        const MatrixX<complex> overlap_ba = overlap_ab.adjoint();
        const MatrixX<complex> overlap_bb = C_beta.adjoint() * S_AO * C_beta;

        // Create the orbital space for the GHF wavefunction model
        const GQCP::OrbitalSpace orbital_space = this->orbitalSpace(this->numberOfSpinors(), this->numberOfElectrons());

        // A KISS implementation of the expectation value of S, from knowdes. (https://gqcg-res.github.io/knowdes/spin-expectation-values-for-ghf.html)
        complex s_x {0.0};
        for (const auto& I : orbital_space.indices(OccupationType::k_occupied)) {  // loop over occupied spinors
            s_x += overlap_ab(I, I).real();
        }

        complex s_y {0.0};
        for (const auto& I : orbital_space.indices(OccupationType::k_occupied)) {  // loop over occupied spinors
            s_y += overlap_ba(I, I).imag();
        }

        complex s_z {0.0};
        for (const auto& I : orbital_space.indices(OccupationType::k_occupied)) {  // loop over occupied spinors
            s_z += 0.5 * (overlap_aa(I, I) - overlap_bb(I, I));
        }

        Vector<complex, 3> s_expectation_value = Vector<complex, 3>::Zero();
        s_expectation_value << s_x, s_y, s_z;

        return s_expectation_value;
    }


    /**
     *  @return the 1-DM expressed in an orthonormal spinor basis related to these optimal GHF parameters
     */
    G1DM<Scalar> calculateOrthonormalBasis1DM() const {

        const auto M = this->numberOfSpinors();
        const auto N = this->numberOfElectrons();
        return GHF<Scalar>::calculateOrthonormalBasis1DM(M, N);
    }


    /**
     *  Construct the partial stability matrix `A` from the GHF stability conditions.
     * 
     *  @note The formula for the `A` matrix is as follows:
     *      A_IAJB = \delta_IJ * F_BA - \delta_BA * F_IJ + (AI||JB)
     * 
     *  @param gsq_hamiltonian      The generalised, second quantized hamiltonian, which contains the necessary two electron operators.
     */
    const GQCP::MatrixX<Scalar> calculatePartialStabilityMatrixA(const GSQHamiltonian<Scalar>& gsq_hamiltonian) const {

        // Create the orbital space to determine the loops.
        const auto orbital_space = this->orbitalSpace();

        // Determine the number of occupied and virtual orbitals.
        const auto& n_occ = orbital_space.numberOfOrbitals(OccupationType::k_occupied);
        const auto& n_virt = orbital_space.numberOfOrbitals(OccupationType::k_virtual);

        // We need the two-electron integrals in MO basis, hence why we transform them with the coefficient matrix.
        // The ground state coefficient matrix is obtained from the QCModel.
        // We need the anti-symmetrized tensor: (AI||JB) = (AI|JB) - (AB|JI). This is obtained by the `.antisymmetrized()` method.
        const auto& g = gsq_hamiltonian.twoElectron().transformed(this->coefficientMatrix()).antisymmetrized();

        // The elements F_BA and F_IJ are the eigenvalues of the one-electron Fock operator. This can be done using the MO energies.
        // The MO energies are contained within the QC Structure.
        std::vector<Scalar> mo_energies;  // We use a std::vector in order to be able to slice the vector later on.

        for (int i = 0; i < this->numberOfSpinors(); i++) {
            mo_energies.push_back(this->orbitalEnergy(i));
        }

        // The MO energies should be split in the respective real and virtual parts.
        std::vector<Scalar> mo_energies_occupied;
        std::copy(mo_energies.begin(), mo_energies.begin() + number_of_occupied_orbitals, std::back_inserter(mo_energies_occupied));

        std::vector<Scalar> mo_energies_virtual;
        std::copy(mo_energies.begin() + number_of_occupied_orbitals, mo_energies.end(), std::back_inserter(mo_energies_virtual));

        // We create a matrix containing the correct elements needed in the formula.
        // These values will be added to the final tensor in the last step of this calculation.
        GQCP::MatrixX<Scalar> F_values(number_of_virtual_orbitals, number_of_occupied_orbitals);

        for (int a = 0; a < number_of_virtual_orbitals; a++) {
            for (int i = 0; i < number_of_occupied_orbitals; i++) {
                F_values(a, i) = mo_energies_virtual[a] + (-1 * mo_energies_occupied[i]);
            }
        }

        // The next step is to create the needed tensor slice.
        GQCP::Tensor<Scalar, 4> A_iajb(number_of_occupied_orbitals, number_of_virtual_orbitals, number_of_occupied_orbitals, number_of_virtual_orbitals);
        for (const auto& i : orbital_space.indices(OccupationType::k_occupied)) {
            for (const auto& a : orbital_space.indices(OccupationType::k_virtual)) {
                for (const auto& j : orbital_space.indices(OccupationType::k_occupied)) {
                    for (const auto& b : orbital_space.indices(OccupationType::k_virtual)) {
                        A_iajb(i, a - number_of_occupied_orbitals, j, b - number_of_occupied_orbitals) = g.parameters()(a, i, j, b);
                    }
                }
            }
        }

        // Add the previously calculated F elements on the correct positions.
        for (int a = 0; a < number_of_virtual_orbitals; a++) {
            for (int i = 0; i < number_of_occupied_orbitals; i++) {
                A_iajb(i, a, i, a) += F_values(a, i);
            }
        }

        // Finally, reshape the tensor to a matrix.
        const GQCP::MatrixX<Scalar> A_matrix = A_iajb.reshape(number_of_occupied_orbitals * number_of_virtual_orbitals, number_of_occupied_orbitals * number_of_virtual_orbitals);

        return A_matrix;
    }


    /**
     *  Construct the partial stability matrix `B` from the GHF stability conditions.
     *
     *  @note The formula for the `B` matrix is as follows:
     *      B_IAJB = (AI||BJ)
     *
     *  @param gsq_hamiltonian      The generalised, second quantized hamiltonian, which contains the necessary two electron operators.
     */
    const GQCP::MatrixX<Scalar> calculatePartialStabilityMatrixB(const GSQHamiltonian<Scalar>& gsq_hamiltonian) const {

        // Create the orbital space to determine the loops.
        const auto orbital_space = this->orbitalSpace();

        // Determine the number of occupied and virtual orbitals.
        const auto& number_of_occupied_orbitals = orbital_space.numberOfOrbitals(OccupationType::k_occupied);
        const auto& number_of_virtual_orbitals = orbital_space.numberOfOrbitals(OccupationType::k_virtual);

        // We need the two-electron integrals in MO basis, hence why we transform them with the coefficient matrix.
        // The ground state coefficient matrix is obtained from the QCModel.
        // We need the anti-symmetrized tensor: (AI||BJ) = (AI|BJ) - (AJ|BI). This is obtained by the `.antisymmetrized()` method.
        const auto& g = gsq_hamiltonian.twoElectron().transformed(this->coefficientMatrix()).antisymmetrized();

        // The next step is to create the needed tensor slice.
        GQCP::Tensor<Scalar, 4> B_iajb(number_of_occupied_orbitals, number_of_virtual_orbitals, number_of_occupied_orbitals, number_of_virtual_orbitals);
        for (const auto& i : orbital_space.indices(OccupationType::k_occupied)) {
            for (const auto& a : orbital_space.indices(OccupationType::k_virtual)) {
                for (const auto& j : orbital_space.indices(OccupationType::k_occupied)) {
                    for (const auto& b : orbital_space.indices(OccupationType::k_virtual)) {
                        B_iajb(i, a - number_of_occupied_orbitals, j, b - number_of_occupied_orbitals) = g.parameters()(a, i, b, j);
                    }
                }
            }
        }

        // Finally, reshape the tensor to a matrix.
        const GQCP::MatrixX<Scalar> B_matrix = B_iajb.reshape(number_of_occupied_orbitals * number_of_virtual_orbitals, number_of_occupied_orbitals * number_of_virtual_orbitals);

        return B_matrix;
    }


    /**
     *  @return the GHF 1-DM in the scalar/AO basis related to these optimal GHF parameters
     */
    G1DM<Scalar> calculateScalarBasis1DM() const {

        const auto N = this->numberOfElectrons();
        return GHF<Scalar>::calculateScalarBasis1DM(this->coefficientMatrix(), N);
    }


    /**
     *  Calculate the GHF stability matrices and return them.
     *
     *  @return The GHF stability matrices.
     */
    GHFStabilityMatrices<Scalar> calculateStabilityMatrices(const GSQHamiltonian<Scalar>& gsq_hamiltonian) const {
        return GHFStabilityMatrices<Scalar> {this->calculatePartialStabilityMatrixA(gsq_hamiltonian), this->calculatePartialStabilityMatrixB(gsq_hamiltonian)};
    }


    /**
     *  @return the coefficient matrix that expresses every spinor orbital (as a column) in the underlying scalar bases
     */
    const GTransformationMatrix<Scalar>& coefficientMatrix() const { return this->C; }

    /**
     *  @return the number of electrons that these GHF model parameters describe
     */
    size_t numberOfElectrons() const { return this->N; }

    /**
     *  @return the number of spinors that these GHF model parameters describe
     */
    size_t numberOfSpinors() const { return this->orbital_energies.size(); }

    /**
     *  @return the orbital energies
     */
    const VectorX<double>& orbitalEnergies() const { return this->orbital_energies; }

    /**
     *  @param i                the index of the orbital
     * 
     *  @return the i-th orbital energy
     */
    double orbitalEnergy(const size_t i) const { return this->orbital_energies(i); }

    /**
     *  @return the implicit occupied-virtual orbital space that is associated to these GHF model parameters
     */
    OrbitalSpace orbitalSpace() const { return GHF<Scalar>::orbitalSpace(this->numberOfSpinors(), this->numberOfElectrons()); }
};

}  // namespace QCModel
}  // namespace GQCP
