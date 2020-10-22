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

        // Convert the matrices Z and D to an Eigen::Tensor<double, 2> D_tensor, as contractions are only implemented for Tensors.
        Eigen::TensorMap<Eigen::Tensor<const Scalar, 2>> P_tensor {P.data(), P.rows(), P.cols()};
        Eigen::TensorMap<Eigen::Tensor<const Scalar, 2>> Z_tensor {Z.parameters().data(), P.rows(), P.cols()};

        // Specify the contraction pair.
        // To calculate the electronic energy, we must perform a double contraction.
        //      0.5 D(mu nu) Z(mu nu)
        // See knowdes: https://gqcg-res.github.io/knowdes/general-hartree-fock-theory.html
        Eigen::array<Eigen::IndexPair<int>, 2> contraction_pair = {Eigen::IndexPair<int>(0, 0), Eigen::IndexPair<int>(1, 1)};

        // Calculate the double contraction (with prefactor 0.5).
        Tensor<Scalar, 0> contraction = 0.5 * P_tensor.contract(Z_tensor, contraction_pair);

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

        Eigen::TensorMap<Eigen::Tensor<const Scalar, 2>> P_aa_tensor {P_aa.data(), P_aa.rows(), P_aa.cols()};
        Eigen::TensorMap<Eigen::Tensor<const Scalar, 2>> P_bb_tensor {P_bb.data(), P_bb.rows(), P_bb.cols()};


        // Specify the contraction pairs for the direct contractions:
        //      P(rho lambda) (mu nu|rho lambda)
        // See knowdes: https://gqcg-res.github.io/knowdes/derivation-of-the-ghf-scf-equations-through-lagrange-multipliers.html
        Eigen::array<Eigen::IndexPair<int>, 2> direct_contraction_pair = {Eigen::IndexPair<int>(0, 2), Eigen::IndexPair<int>(1, 3)};

        // Do the actual contractions, and convert the given tensor back to a matrix.
        const auto& g = sq_hamiltonian.twoElectron().parameters();
        Tensor<Scalar, 2> J_tensor = P_aa_tensor.contract(g.Eigen(), direct_contraction_pair) + P_bb_tensor.contract(g.Eigen(), direct_contraction_pair);

        Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>> J {J_tensor.data(), J_tensor.dimension(0), J_tensor.dimension(1)};

        return ScalarGSQOneElectronOperator<Scalar>(J);
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


        Eigen::TensorMap<Eigen::Tensor<const Scalar, 2>> P_aa_tensor {P_aa.data(), P_aa.rows(), P_aa.cols()};
        Eigen::TensorMap<Eigen::Tensor<const Scalar, 2>> P_ab_tensor {P_ab.data(), P_ab.rows(), P_ab.cols()};
        Eigen::TensorMap<Eigen::Tensor<const Scalar, 2>> P_ba_tensor {P_ba.data(), P_ba.rows(), P_ba.cols()};
        Eigen::TensorMap<Eigen::Tensor<const Scalar, 2>> P_bb_tensor {P_bb.data(), P_bb.rows(), P_bb.cols()};


        // Specify the contraction pairs for the exchange contractions:
        //      P(lambda rho) (mu rho|lambda nu)
        Eigen::array<Eigen::IndexPair<int>, 2> exchange_contraction_pair = {Eigen::IndexPair<int>(0, 2), Eigen::IndexPair<int>(1, 1)};

        // Do the actual contractions, and convert the given tensor back to a matrix.
        const auto& g = sq_hamiltonian.twoElectron().parameters();

        Tensor<Scalar, 2> K_aa_tensor = P_aa_tensor.contract(g.Eigen(), exchange_contraction_pair);
        Tensor<Scalar, 2> K_ab_tensor = P_ba_tensor.contract(g.Eigen(), exchange_contraction_pair);
        Tensor<Scalar, 2> K_ba_tensor = P_ab_tensor.contract(g.Eigen(), exchange_contraction_pair);
        Tensor<Scalar, 2> K_bb_tensor = P_bb_tensor.contract(g.Eigen(), exchange_contraction_pair);

        Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>> K_aa {K_aa_tensor.data(), K_aa_tensor.dimension(0), K_aa_tensor.dimension(1)};
        Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>> K_ab {K_ab_tensor.data(), K_ab_tensor.dimension(0), K_ab_tensor.dimension(1)};
        Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>> K_ba {K_ba_tensor.data(), K_ba_tensor.dimension(0), K_ba_tensor.dimension(1)};
        Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>> K_bb {K_bb_tensor.data(), K_bb_tensor.dimension(0), K_bb_tensor.dimension(1)};


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
     *  @return the GHF 1-DM in the scalar/AO basis related to these optimal GHF parameters
     */
    G1DM<Scalar> calculateScalarBasis1DM() const {

        const auto N = this->numberOfElectrons();
        return GHF<Scalar>::calculateScalarBasis1DM(this->coefficientMatrix(), N);
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
};  // namespace QCModel


}  // namespace QCModel
}  // namespace GQCP
