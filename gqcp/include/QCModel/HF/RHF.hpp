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


#include "Basis/SpinorBasis/OrbitalSpace.hpp"
#include "Basis/Transformations/RTransformation.hpp"
#include "DensityMatrix/Orbital1DM.hpp"
#include "Mathematical/Representation/ImplicitRankFourTensorSlice.hpp"
#include "Mathematical/Representation/SquareMatrix.hpp"
#include "Operator/SecondQuantized/RSQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCModel/HF/StabilityMatrices/RHFStabilityMatrices.hpp"
#include "QuantumChemical/Spin.hpp"


namespace GQCP {
namespace QCModel {


/**
 *  The restricted Hartree-Fock wave function model.
 * 
 *  @tparam _Scalar             the type of scalar that is used for the expansion of the spatial orbitals in their underlying scalar basis
 */
template <typename _Scalar>
class RHF {
public:
    using Scalar = _Scalar;


private:
    size_t N_P;  // the number of electron pairs

    VectorX<double> orbital_energies;  // sorted in ascending energies

    // The transformation that expresses the RHF MOs in terms of the atomic spinors.
    RTransformation<Scalar> C;


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  The standard member-wise constructor
     * 
     *  @param N_P                  the number of electron pairs
     *  @param C                    The transformation that expresses the RHF MOs in terms of the atomic spinors.
     *  @param orbital_energies     the RHF MO energies
     */
    RHF(const size_t N_P, const VectorX<double>& orbital_energies, const RTransformation<double>& C) :
        N_P {N_P},
        orbital_energies {orbital_energies},
        C(C) {}


    /**
     *  Default constructor setting everything to zero
     */
    RHF() :
        RHF(0.0, RTransformation<double>::Zero(0), VectorX<double>::Zero(0)) {}


    /*
     *  STATIC PUBLIC METHODS
     */


    /**
     *  @param D                the RHF density matrix in a scalar basis
     *  @param H_core           the core Hamiltonian expressed in the same scalar basis
     *  @param F                the Fock matrix in the same scalar basis
     *
     *  @return the RHF electronic energy
     */
    static double calculateElectronicEnergy(const Orbital1DM<Scalar>& D, const ScalarRSQOneElectronOperator<Scalar>& H_core, const ScalarRSQOneElectronOperator<Scalar>& F) {

        // First, calculate the sum of H_core and F (this saves a contraction).
        const auto Z = H_core + F;

        // Convert the matrix Z to an GQCP::Tensor<double, 2> Z_tensor.
        // Einsum is only implemented for a tensor + a matrix, not for 2 matrices.
        Eigen::TensorMap<Eigen::Tensor<const Scalar, 2>> Z_t {Z.parameters().data(), D.rows(), D.cols()};
        Tensor<Scalar, 2> Z_tensor = Tensor<Scalar, 2>(Z_t);

        // To calculate the electronic energy, we must perform a double contraction (with prefactor 0.5)
        //      0.5 D(nu mu) Z(mu nu)
        Tensor<Scalar, 0> contraction = 0.5 * Z_tensor.template einsum<2>("ij,ji->", D);

        // As the double contraction of two matrices is a scalar (a tensor of rank 0), we should access the value as (0)
        return contraction(0);
    }


    /**
     *  @param F                the Fock matrix expressed in a scalar basis
     *  @param D                the RHF density matrix in the same scalar basis
     *  @param S                the overlap matrix of that scalar basis
     * 
     *  @return the RHF error matrix
     */
    static SquareMatrix<Scalar> calculateError(const SquareMatrix<Scalar>& F, const Orbital1DM<Scalar>& D, const SquareMatrix<Scalar>& S) {
        return F * D * S - S * D * F;
    }


    /**
     *  @param sq_hamiltonian       the Hamiltonian expressed in an orthonormal basis
     *  @param N_P                  the number of electron pairs
     *  @param a                    the first virtual orbital index
     *  @param i                    the first occupied orbital index
     *  @param b                    the second virtual orbital index
     *  @param j                    the second occupied orbital index
     * 
     *  @return an element of the RHF orbital Hessian
     */
    static Scalar calculateOrbitalHessianElement(const RSQHamiltonian<Scalar>& sq_hamiltonian, const size_t N_P, const size_t a, const size_t i, const size_t b, const size_t j) {

        // Prepare some variables.
        const auto& g_op = sq_hamiltonian.twoElectron();
        const auto K = g_op.numberOfOrbitals();  // The number of spatial orbitals.
        const auto& g = g_op.parameters();

        const auto orbital_space = RHF<Scalar>::orbitalSpace(K, N_P);


        double value {0.0};

        // Inactive Fock matrix part
        const auto F = sq_hamiltonian.calculateInactiveFockian(orbital_space).parameters();
        if (i == j) {
            value += F(a, b);
        }

        if (a == b) {
            value -= F(i, j);
        }


        // Two-electron part
        value += 4 * g(a, i, b, j) - g(a, b, i, j) - g(a, j, b, i);

        return 4 * value;
    }


    /**
     *  @param sq_hamiltonian       the Hamiltonian expressed in an orthonormal basis
     *  @param N_P                  the number of electron pairs
     * 
     *  @return the RHF orbital Hessian as a ImplicitRankFourTensorSlice, i.e. an object with a suitable operator() implemented
     */
    static ImplicitRankFourTensorSlice<Scalar> calculateOrbitalHessianTensor(const RSQHamiltonian<Scalar>& sq_hamiltonian, const size_t N_P) {

        // Create an occupied-virtual orbital space.
        const auto K = sq_hamiltonian.numberOfOrbitals();
        const auto orbital_space = OrbitalSpace::Implicit({{OccupationType::k_occupied, N_P}, {OccupationType::k_virtual, K - N_P}});  // N_P occupied spatial orbitals, K-N_P virtual spatial orbitals

        auto hessian = orbital_space.template initializeRepresentableObjectFor<Scalar>(OccupationType::k_virtual, OccupationType::k_occupied, OccupationType::k_virtual, OccupationType::k_occupied);  // zero-initialize an object suitable for the representation of a virtual-occupied,virtual-occupied object (ai,bj)

        // Loop over all indices (ai,bj) to construct the orbital hessian
        for (const auto& a : orbital_space.indices(OccupationType::k_virtual)) {
            for (const auto& i : orbital_space.indices(OccupationType::k_occupied)) {

                for (const auto& b : orbital_space.indices(OccupationType::k_virtual)) {
                    for (const auto& j : orbital_space.indices(OccupationType::k_occupied)) {
                        hessian(a, i, b, j) = RHF<Scalar>::calculateOrbitalHessianElement(sq_hamiltonian, N_P, a, i, b, j);
                    }
                }
            }
        }

        return hessian;
    }


    /**
     *  Calculate the RHF 1-DM expressed in an orthonormal spin-orbital basis.
     * 
     *  @param K            The number of spatial orbitals.
     *  @param N            The total number of electrons.
     *
     *  @return The RHF 1-DM expressed in an orthonormal spin-orbital basis.
     */
    static Orbital1DM<Scalar> calculateOrthonormalBasis1DM(const size_t K, const size_t N) {

        if (N % 2 != 0) {
            throw std::invalid_argument("QCMethod::RHF::calculateOrthonormalBasis1DM(const size_t, const size_t): The number of given electrons cannot be odd for RHF.");
        }

        // The 1-DM for RHF looks like (for K=5, N=6)
        //    2  0  0  0  0
        //    0  2  0  0  0
        //    0  0  2  0  0
        //    0  0  0  0  0
        //    0  0  0  0  0

        Orbital1DM<double> D = Orbital1DM<double>::Zero(K);
        D.topLeftCorner(N / 2, N / 2) = 2 * SquareMatrix<double>::Identity(N / 2);

        return D;
    }


    /**
     *  Calculate the RHF 2-DM expressed in an orthonormal spin-orbital basis.
     * 
     *  @param K            The number of spatial orbitals.
     *  @param N            The total number of electrons.
     *
     *  @return The RHF 2-DM expressed in an orthonormal spin-orbital basis.
     */
    static Orbital2DM<Scalar> calculateOrthonormalBasis2DM(const size_t K, const size_t N) {

        if (N % 2 != 0) {
            throw std::invalid_argument("QCMethod::RHF::calculateOrthonormalBasis1DM(const size_t, const size_t): The number of given electrons cannot be odd for RHF.");
        }

        const size_t N_P = N / 2;  // The number of electron pairs.


        // Create the orbital space to determine the loops.
        const auto orbital_space = RHF<Scalar>::orbitalSpace(K, N_P);


        // Implement a KISS formula for the RHF 2-DM.
        Orbital2DM<Scalar> d = Orbital2DM<Scalar>::Zero(K);

        for (const auto& i : orbital_space.indices(OccupationType::k_occupied)) {
            for (const auto& j : orbital_space.indices(OccupationType::k_occupied)) {
                for (const auto& k : orbital_space.indices(OccupationType::k_occupied)) {
                    for (const auto& l : orbital_space.indices(OccupationType::k_occupied)) {
                        if ((i == j) && (k == l)) {
                            d(i, j, k, l) += 4.0;
                        }

                        if ((i == l) && (j == k)) {
                            d(i, j, k, l) -= 2.0;
                        }
                    }
                }
            }
        }

        return d;
    }


    /**
     *  @param C    the coefficient matrix that expresses every spatial orbital (as a column) in its underlying scalar basis
     *  @param N    the number of electrons
     *
     *  @return the RHF 1-DM expressed in the underlying scalar basis
     */
    static Orbital1DM<Scalar> calculateScalarBasis1DM(const RTransformation<double>& C, const size_t N) {

        const size_t K = C.numberOfOrbitals();
        const auto D_orthonormal = RHF<Scalar>::calculateOrthonormalBasis1DM(K, N);

        // Transform the 1-DM in an orthonormal basis to the underlying scalar basis.
        return D_orthonormal.transformed(C.inverse());
    }


    /**
     *  Calculate the RHF Fock matrix F = H_core + G, in which G is a contraction of the density matrix and the two-electron integrals
     *
     *  @param D                    the RHF density matrix in a scalar basis
     *  @param sq_hamiltonian       the Hamiltonian expressed in the same scalar basis
     *
     *  @return the RHF Fock matrix expressed in the scalar basis
     */
    static ScalarRSQOneElectronOperator<Scalar> calculateScalarBasisFockMatrix(const Orbital1DM<Scalar>& D, const RSQHamiltonian<Scalar>& sq_hamiltonian) {

        // get the two-electron parameters
        const auto& g = sq_hamiltonian.twoElectron().parameters();

        // To calculate G, we must perform two double contractions
        //      1. (mu nu|rho lambda) P(lambda rho)
        const Tensor<Scalar, 2> direct_contraction = g.template einsum<2>("ijkl,lk->ij", D);
        //      2. -0.5 (mu lambda|rho nu) P(lambda rho)
        const Tensor<Scalar, 2> exchange_contraction = -0.5 * g.template einsum<2>("ilkj,lk->ij", D);

        // The previous contractions are Tensor<Scalar, 2> instances. In order to calculate the total G matrix, we will convert them back into GQCP::Matrix<Scalar>
        auto G1 = direct_contraction.asMatrix();
        auto G2 = exchange_contraction.asMatrix();

        return ScalarRSQOneElectronOperator<Scalar> {sq_hamiltonian.core().parameters() + G1 + G2};
    }


    /**
     *  @param N            the number of electrons
     *
     *  @return the (spatial orbital, not spin-orbital) index of the RHF HOMO in an implicit orbital space
     */
    static size_t homoIndex(const size_t N) {

        if (N % 2 != 0) {
            throw std::invalid_argument("QCModel::RHF::homoIndex(const size_t): Can't calculate the RHF HOMO index for an odd number of electrons N.");
        }

        return N / 2 - 1;  // need to subtract 1 because computer indices start at 0
    }


    /**
     *  @param K            the number of spatial orbitals
     *  @param N            the number of electrons
     *
     *  @return the (spatial orbital, not spin-orbital) index of the RHF LUMO in an implicit orbital space
     */
    static size_t lumoIndex(const size_t K, const size_t N) {

        if (N >= 2 * K) {
            throw std::invalid_argument("QCModel::RHF::lumoIndex(const size_t, constsize_t): There is no LUMO for the given number of electrons N and spatial orbitals K");
        }

        return RHF<Scalar>::homoIndex(N) + 1;
    }


    /**
     *  @param K            the number of spatial orbitals
     *  @param N_P          the number of electrons
     * 
     *  @return the implicit (i.e. with ascending and contiguous orbital indices) occupied-virtual orbital space that corresponds to these RHF model parameters
     */
    static OrbitalSpace orbitalSpace(const size_t K, const size_t N_P) {

        return OrbitalSpace::Implicit({{OccupationType::k_occupied, N_P}, {OccupationType::k_virtual, K - N_P}});
    }


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return The 1-DM expressed in an orthonormal spinor basis related to these optimal RHF parameters.
     */
    Orbital1DM<Scalar> calculateOrthonormalBasis1DM() const {

        const auto K = this->numberOfSpatialOrbitals();
        const auto N = 2 * this->numberOfElectronPairs();
        return RHF<Scalar>::calculateOrthonormalBasis1DM(K, N);
    }


    /**
     *  @return The 2-DM expressed in an orthonormal spinor basis related to these optimal RHF parameters.
     */
    Orbital2DM<Scalar> calculateOrthonormalBasis2DM() const {

        const auto K = this->numberOfSpatialOrbitals();
        const auto N = 2 * this->numberOfElectronPairs();
        return RHF<Scalar>::calculateOrthonormalBasis2DM(K, N);
    }


    /**
     *  @return A matrix containing all the possible excitation energies of the wavefunction model.
     * 
     *  @note The rows are determined by the number of virtual orbitals, the columns by the number of occupied orbitals.
     */
    GQCP::MatrixX<Scalar> excitationEnergies() const {

        // Create the orbital space to determine the loops.
        const auto orbital_space = this->orbitalSpace();

        // Determine the number of occupied and virtual orbitals.
        const auto& n_occ = orbital_space.numberOfOrbitals(OccupationType::k_occupied);
        const auto& n_virt = orbital_space.numberOfOrbitals(OccupationType::k_virtual);

        // Calculate the occupied and virtual orbital energies.
        const auto occupied_energies = this->occupiedOrbitalEnergies();
        const auto virtual_energies = this->virtualOrbitalEnergies();

        // Create the F matrix.
        GQCP::MatrixX<Scalar> F_values(n_virt, n_occ);
        for (int a = 0; a < n_virt; a++) {
            for (int i = 0; i < n_occ; i++) {
                F_values(a, i) = virtual_energies[a] - occupied_energies[i];
            }
        }
        return F_values;
    }


    /**
     *  @return the RHF 1-DM in the scalar/AO basis related to these optimal RHF parameters
     */
    Orbital1DM<Scalar> calculateScalarBasis1DM() const {

        const auto N = 2 * this->numberOfElectronPairs();
        return RHF<Scalar>::calculateScalarBasis1DM(this->expansion(), N);
    }


    /**
     *  Construct the `singlet A` stability matrix from the RHF stability conditions.
     * 
     *  @note The formula for the `singlet A` matrix is as follows:
     *      A_IAJB = \delta_IJ * (F_R)_BA - \delta_AB * (F_R)_IJ + 2 * (AI|JB) - (AB|JI).
     * 
     *  @param rsq_hamiltonian      The second quantized Hamiltonian, expressed in the orthonormal, 'restricted' spin orbital basis of the RHF MOs, which contains the necessary two-electron operators.
     * 
     *  @return The singlet-A stability matrix.
     */
    GQCP::MatrixX<Scalar> calculateSingletAStabilityMatrix(const RSQHamiltonian<Scalar>& rsq_hamiltonian) const {

        // Create the orbital space.
        const auto orbital_space = this->orbitalSpace();

        // Create the number of occupied and virtual orbitals.
        const auto& n_occ = orbital_space.numberOfOrbitals(OccupationType::k_occupied);
        const auto& n_virt = orbital_space.numberOfOrbitals(OccupationType::k_virtual);

        // The two electron integrals are extracted from the Hamiltonian.
        const auto& g = rsq_hamiltonian.twoElectron().parameters();

        // The elements (F_R)_AA and (F_R)_IJ are the eigenvalues of the one-electron Fock operator.
        // The excitationEnergies API can be used to find these values.
        const auto F_values = this->excitationEnergies();

        // The next step is to create the needed tensor slice.
        // Zero-initialize an occupied-virtual-occupied-virtual object.
        auto singlet_A_slice = orbital_space.template initializeRepresentableObjectFor<Scalar>(OccupationType::k_occupied, OccupationType::k_virtual, OccupationType::k_occupied, OccupationType::k_virtual);
        for (const auto& i : orbital_space.indices(OccupationType::k_occupied)) {
            for (const auto& a : orbital_space.indices(OccupationType::k_virtual)) {
                for (const auto& j : orbital_space.indices(OccupationType::k_occupied)) {
                    for (const auto& b : orbital_space.indices(OccupationType::k_virtual)) {
                        singlet_A_slice(i, a, j, b) = 2 * g(a, i, j, b) - g(a, b, j, i);
                    }
                }
            }
        }

        // Turn the ImplicitRankFourTensorSlice in an actual Tensor.
        auto singlet_A_iajb = singlet_A_slice.asTensor();

        // Add the previously calculated F values on the correct positions.
        for (int a = 0; a < n_virt; a++) {
            for (int i = 0; i < n_occ; i++) {
                singlet_A_iajb(i, a, i, a) += F_values(a, i);
            }
        }

        // Finally, reshape the tensor to a matrix.
        const auto singlet_A_matrix = singlet_A_iajb.reshape(n_occ * n_virt, n_occ * n_virt);

        return singlet_A_matrix;
    }


    /**
     *  Construct the `singlet B` stability matrix from the RHF stability conditions.
     * 
     *  @note The formula for the `singlet B` matrix is as follows:
     *      B_IAJB = 2 * (AI|BJ) - (AJ|BI).
     * 
     *  @param rsq_hamiltonian      The second quantized Hamiltonian, expressed in the orthonormal, 'restricted' spin orbital basis of the RHF MOs, which contains the necessary two-electron operators.
     * 
     *  @return The singlet-B stability matrix.
     */
    GQCP::MatrixX<Scalar> calculateSingletBStabilityMatrix(const RSQHamiltonian<Scalar>& rsq_hamiltonian) const {

        // Create the orbital space.
        const auto orbital_space = this->orbitalSpace();

        // Create the number of occupied and virtual orbitals.
        const auto& n_occ = orbital_space.numberOfOrbitals(OccupationType::k_occupied);
        const auto& n_virt = orbital_space.numberOfOrbitals(OccupationType::k_virtual);

        // The two electron integrals are extracted from the Hamiltonian.
        const auto& g = rsq_hamiltonian.twoElectron().parameters();

        // The next step is to create the needed tensor slice.
        // Zero-initialize an occupied-virtual-occupied-virtual object.
        auto singlet_B_slice = orbital_space.template initializeRepresentableObjectFor<Scalar>(OccupationType::k_occupied, OccupationType::k_virtual, OccupationType::k_occupied, OccupationType::k_virtual);
        for (const auto& i : orbital_space.indices(OccupationType::k_occupied)) {
            for (const auto& a : orbital_space.indices(OccupationType::k_virtual)) {
                for (const auto& j : orbital_space.indices(OccupationType::k_occupied)) {
                    for (const auto& b : orbital_space.indices(OccupationType::k_virtual)) {
                        singlet_B_slice(i, a, j, b) = 2 * g(a, i, b, j) - g(a, j, b, i);
                    }
                }
            }
        }

        // Turn the ImplicitRankFourTensorSlice in an actual Tensor.
        auto singlet_B_iajb = singlet_B_slice.asTensor();

        // Finally, reshape the tensor to a matrix.
        const auto singlet_B_matrix = singlet_B_iajb.reshape(n_occ * n_virt, n_occ * n_virt);

        return singlet_B_matrix;
    }


    /**
     *  Construct the `triplet A` stability matrix from the RHF stability conditions.
     * 
     *  @note The formula for the `triplet A` matrix is as follows:
     *      A_IAJB = \delta_IJ * (F_R)_BA - \delta_AB * (F_R)_IJ - (AB|JI).
     * 
     *  @param rsq_hamiltonian      The second quantized Hamiltonian, expressed in the orthonormal, 'restricted' spin orbital basis of the RHF MOs, which contains the necessary two-electron operators.
     * 
     *  @return the triplet-A stability matrix.
     */
    GQCP::MatrixX<Scalar> calculateTripletAStabilityMatrix(const RSQHamiltonian<Scalar>& rsq_hamiltonian) const {

        // Create the orbital space.
        const auto orbital_space = this->orbitalSpace();

        // Create the number of occupied and virtual orbitals.
        const auto& n_occ = orbital_space.numberOfOrbitals(OccupationType::k_occupied);
        const auto& n_virt = orbital_space.numberOfOrbitals(OccupationType::k_virtual);

        // The two electron integrals are extracted from the Hamiltonian.
        const auto& g = rsq_hamiltonian.twoElectron().parameters();

        // The elements (F_R)_AA and (F_R)_IJ are the eigenvalues of the one-electron Fock operator.
        // The excitationEnergies API can be used to find these values.
        const auto F_values = this->excitationEnergies();

        // The next step is to create the needed tensor slice.
        // Zero-initialize an occupied-virtual-occupied-virtual object.
        auto triplet_A_slice = orbital_space.template initializeRepresentableObjectFor<Scalar>(OccupationType::k_occupied, OccupationType::k_virtual, OccupationType::k_occupied, OccupationType::k_virtual);
        for (const auto& i : orbital_space.indices(OccupationType::k_occupied)) {
            for (const auto& a : orbital_space.indices(OccupationType::k_virtual)) {
                for (const auto& j : orbital_space.indices(OccupationType::k_occupied)) {
                    for (const auto& b : orbital_space.indices(OccupationType::k_virtual)) {
                        triplet_A_slice(i, a, j, b) = -g(a, b, j, i);
                    }
                }
            }
        }

        // Turn the ImplicitRankFourTensorSlice in an actual Tensor.
        auto triplet_A_iajb = triplet_A_slice.asTensor();

        // Add the previously calculated F values on the correct positions.
        for (int a = 0; a < n_virt; a++) {
            for (int i = 0; i < n_occ; i++) {
                triplet_A_iajb(i, a, i, a) += F_values(a, i);
            }
        }

        // Finally, reshape the tensor to a matrix.
        const auto triplet_A_matrix = triplet_A_iajb.reshape(n_occ * n_virt, n_occ * n_virt);

        return triplet_A_matrix;
    }


    /**
     *  Construct the `triplet B` stability matrix from the RHF stability conditions.
     * 
     *  @note The formula for the `triplet B` matrix is as follows:
     *      B_IAJB = - (AJ|BI).
     * 
     *  @param rsq_hamiltonian      The second quantized Hamiltonian, expressed in the orthonormal, 'restricted' spin orbital basis of the RHF MOs, which contains the necessary two-electron operators.
     * 
     *  @return The triplet-B stability matrix.
     */
    GQCP::MatrixX<Scalar> calculateTripletBStabilityMatrix(const RSQHamiltonian<Scalar>& rsq_hamiltonian) const {

        // Create the orbital space.
        const auto orbital_space = this->orbitalSpace();

        // Create the number of occupied and virtual orbitals.
        const auto& n_occ = orbital_space.numberOfOrbitals(OccupationType::k_occupied);
        const auto& n_virt = orbital_space.numberOfOrbitals(OccupationType::k_virtual);

        // The two electron integrals are extracted from the Hamiltonian.
        const auto& g = rsq_hamiltonian.twoElectron().parameters();

        // The next step is to create the needed tensor slice.
        // Zero-initialize an occupied-virtual-occupied-virtual object.
        auto triplet_B_slice = orbital_space.template initializeRepresentableObjectFor<Scalar>(OccupationType::k_occupied, OccupationType::k_virtual, OccupationType::k_occupied, OccupationType::k_virtual);
        for (const auto& i : orbital_space.indices(OccupationType::k_occupied)) {
            for (const auto& a : orbital_space.indices(OccupationType::k_virtual)) {
                for (const auto& j : orbital_space.indices(OccupationType::k_occupied)) {
                    for (const auto& b : orbital_space.indices(OccupationType::k_virtual)) {
                        triplet_B_slice(i, a, j, b) = -g(a, j, b, i);
                    }
                }
            }
        }

        // Turn the ImplicitRankFourTensorSlice in an actual Tensor.
        auto triplet_B_iajb = triplet_B_slice.asTensor();

        // Finally, reshape the tensor to a matrix.
        const auto triplet_B_matrix = triplet_B_iajb.reshape(n_occ * n_virt, n_occ * n_virt);

        return triplet_B_matrix;
    }


    /**
     *  Calculate the RHF stability matrices and return them.
     * 
     *  @param rsq_hamiltonian      The second quantized Hamiltonian, expressed in the orthonormal, 'restricted' spin orbital basis of the RHF MOs, which contains the necessary two-electron operators.
     *
     *  @return The RHF stability matrices.
     */
    RHFStabilityMatrices<Scalar> calculateStabilityMatrices(const RSQHamiltonian<Scalar>& rsq_hamiltonian) const {
        return RHFStabilityMatrices<Scalar> {this->calculateSingletAStabilityMatrix(rsq_hamiltonian), this->calculateSingletBStabilityMatrix(rsq_hamiltonian),
                                             this->calculateTripletAStabilityMatrix(rsq_hamiltonian), this->calculateTripletBStabilityMatrix(rsq_hamiltonian)};
    }


    /**
     *  @return The transformation that expresses the RHF MOs in terms of the underlying AOs.
     */
    const RTransformation<Scalar>& expansion() const { return this->C; }

    /**
     *  @param N            the number of electrons
     *
     *  @return the (spatial orbital, not spin-orbital) index of the RHF HOMO in an implicit orbital space
     */
    size_t homoIndex() const { return RHF<Scalar>::homoIndex(this->numberOfElectrons()); }

    /**
     *  @param K            the number of spatial orbitals
     *  @param N            the number of electrons
     *
     *  @return the (spatial orbital, not spin-orbital) index of the RHF LUMO in an implicit orbital space
     */
    size_t lumoIndex() const { return RHF<Scalar>::lumoIndex(this->numberOfSpatialOrbitals(), this->numberOfElectrons()); }

    /**
     *  @return the number of electron pairs that these RHF model parameters describe
     */
    size_t numberOfElectronPairs() const { return this->N_P; }

    /**
     *  @return the total number of electrons that these RHF model parameters describe
     */
    size_t numberOfElectrons() const { return 2 * this->numberOfElectronPairs(); }

    /**
     *  @param sigma            the spin of the electrons (alpha or beta)
     * 
     *  @return the number of sigma-electrons that these RHF model parameters describe
     */
    size_t numberOfElectrons(const Spin sigma) const { return this->numberOfElectronPairs(); }

    /**
     *  @return the number of spatial orbitals that these RHF model parameters describe
     */
    size_t numberOfSpatialOrbitals() const { return this->expansion().numberOfOrbitals(); }

    /**
     *  @return The orbital energies belonging to the occupied orbitals.
     */
    std::vector<double> occupiedOrbitalEnergies() const {

        // Determine the number of occupied orbitals.
        const auto n_occ = this->orbitalSpace().numberOfOrbitals(OccupationType::k_occupied);

        std::vector<double> mo_energies;  // We use a std::vector in order to be able to slice the vector later on.
        for (int i = 0; i < this->numberOfSpatialOrbitals(); i++) {
            mo_energies.push_back(this->orbitalEnergy(i));
        }

        // Add the values with indices smaller than the occupied orbital indices, to the new vector.
        std::vector<double> mo_energies_occupied;
        std::copy(mo_energies.begin(), mo_energies.begin() + n_occ, std::back_inserter(mo_energies_occupied));
        return mo_energies_occupied;
    }


    /**
     *  @return all the spatial orbital energies
     */
    const VectorX<double>& orbitalEnergies() const { return this->orbital_energies; }

    /**
     *  @param i                the index of the orbital
     * 
     *  @return the i-th orbital energy
     */
    double orbitalEnergy(const size_t i) const { return this->orbital_energies(i); }

    /**
     *  @return the implicit occupied-virtual orbital space that is associated to these RHF model parameters
     */
    OrbitalSpace orbitalSpace() const { return RHF<Scalar>::orbitalSpace(this->numberOfSpatialOrbitals(), this->numberOfElectronPairs()); }

    /**
     *  @return all the spin-orbital energies, with the alpha spin-orbital energies appearing before the beta spin-orbital energies
     */
    VectorX<double> spinOrbitalEnergiesBlocked() const {

        const auto K = this->numberOfSpatialOrbitals();

        GQCP::VectorX<double> total_orbital_energies {2 * K};
        total_orbital_energies.head(K) = this->orbitalEnergies();
        total_orbital_energies.tail(K) = this->orbitalEnergies();

        return total_orbital_energies;
    }


    /**
     *  @return all the spin-orbital energies, with the alpha and beta-spinorbital energies interleaved
     */
    VectorX<double> spinOrbitalEnergiesInterleaved() const {

        const auto K = this->numberOfSpatialOrbitals();

        GQCP::VectorX<double> total_orbital_energies {2 * K};
        for (size_t p = 0; p < K; p++) {
            total_orbital_energies(2 * p) = this->orbitalEnergy(p);
            total_orbital_energies(2 * p + 1) = this->orbitalEnergy(p);
        }

        return total_orbital_energies;
    }


    /**
     *  @return The orbital energies belonging to the virtual orbitals.
     */
    std::vector<double> virtualOrbitalEnergies() const {

        // Determine the number of occupied orbitals.
        const auto n_occ = this->orbitalSpace().numberOfOrbitals(OccupationType::k_occupied);

        std::vector<double> mo_energies;  // We use a std::vector in order to be able to slice the vector later on.
        for (int i = 0; i < this->numberOfSpatialOrbitals(); i++) {
            mo_energies.push_back(this->orbitalEnergy(i));
        }

        // Add the values with indices greater than the occupied orbital indices, i.e. the virtual orbital indices, to the new vector.
        std::vector<double> mo_energies_virtual;
        std::copy(mo_energies.begin() + n_occ, mo_energies.end(), std::back_inserter(mo_energies_virtual));
        return mo_energies_virtual;
    }

};  // namespace QCModel


}  // namespace QCModel
}  // namespace GQCP
