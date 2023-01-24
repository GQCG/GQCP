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


#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Basis/SpinorBasis/CurrentDensityMatrixElement.hpp"
#include "Basis/SpinorBasis/OrbitalSpace.hpp"
#include "Basis/Transformations/RTransformation.hpp"
#include "DensityMatrix/Orbital1DM.hpp"
#include "Mathematical/Grid/CubicGrid.hpp"
#include "Mathematical/Representation/ImplicitRankFourTensorSlice.hpp"
#include "Mathematical/Representation/LeviCivitaTensor.hpp"
#include "Mathematical/Representation/SquareMatrix.hpp"
#include "Operator/SecondQuantized/EvaluableRSQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/RSQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCModel/HF/StabilityMatrices/RHFStabilityMatrices.hpp"
#include "QuantumChemical/Spin.hpp"
#include "Utilities/complex.hpp"
#include "Utilities/type_traits.hpp"


namespace GQCP {
namespace QCModel {


/**
 *  The restricted Hartree-Fock wave function model.
 *
 *  @tparam _Scalar             The type of scalar that is used for the expansion of the spatial orbitals in their underlying scalar basis: real or complex.
 */
template <typename _Scalar>
class RHF {
public:
    // The type of scalar that is used for the expansion of the spatial orbitals in their underlying scalar basis: real or complex.
    using Scalar = _Scalar;


private:
    // The number of electron pairs.
    size_t N_P;

    // The orbital energies sorted in ascending order.
    VectorX<Scalar> orbital_energies;

    // The transformation that expresses the RHF MOs in terms of the atomic spinors.
    RTransformation<Scalar> C;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  The standard member-wise constructor.
     *
     *  @param N_P                  The number of electron pairs.
     *  @param C                    The transformation that expresses the RHF MOs in terms of the atomic spinors.
     *  @param orbital_energies     The RHF MO energies.
     */
    RHF(const size_t N_P, const VectorX<Scalar>& orbital_energies, const RTransformation<Scalar>& C) :
        N_P {N_P},
        orbital_energies {orbital_energies},
        C(C) {}


    /**
     *  Default constructor setting everything to zero.
     */
    RHF() :
        RHF(0.0, RTransformation<Scalar>::Zero(0), VectorX<Scalar>::Zero(0)) {}


    /*
     *  MARK: Energy
     */

    /**
     *  @param D                The RHF density matrix in a scalar basis.
     *  @param H_core           The core Hamiltonian expressed in the same scalar basis.
     *  @param F                The Fock matrix in the same scalar basis.
     *
     *  @return The RHF electronic energy.
     */
    static Scalar calculateElectronicEnergy(const Orbital1DM<Scalar>& D, const ScalarRSQOneElectronOperator<Scalar>& H_core, const ScalarRSQOneElectronOperator<Scalar>& F) {

        // First, calculate the sum of H_core and F (this saves a contraction).
        const auto Z = H_core + F;

        // Convert the matrix Z to an GQCP::Tensor<double, 2> Z_tensor.
        // Einsum is only implemented for a tensor + a matrix, not for 2 matrices.
        Eigen::TensorMap<Eigen::Tensor<const Scalar, 2>> Z_t {Z.parameters().data(), D.matrix().rows(), D.matrix().cols()};
        Tensor<Scalar, 2> Z_tensor = Tensor<Scalar, 2>(Z_t);

        // To calculate the electronic energy, we must perform a double contraction (with prefactor 0.5):
        //      0.5 D(nu mu) Z(mu nu).
        Tensor<Scalar, 0> contraction = 0.5 * Z_tensor.template einsum<2>("ij,ij->", D.matrix());

        // As the double contraction of two matrices is a scalar (a tensor of rank 0), we should access the value as (0).
        return contraction(0);
    }


    /**
     *  @return A matrix containing all the possible excitation energies of the wavefunction model.
     *
     *  @note The rows are determined by the number of virtual orbitals, the columns by the number of occupied orbitals.
     */
    MatrixX<Scalar> excitationEnergies() const {

        // Create the orbital space to determine the loops.
        const auto orbital_space = this->orbitalSpace();

        // Determine the number of occupied and virtual orbitals.
        const auto& n_occ = orbital_space.numberOfOrbitals(OccupationType::k_occupied);
        const auto& n_virt = orbital_space.numberOfOrbitals(OccupationType::k_virtual);

        // Calculate the occupied and virtual orbital energies.
        const auto occupied_energies = this->occupiedOrbitalEnergies();
        const auto virtual_energies = this->virtualOrbitalEnergies();

        // Create the F matrix.
        MatrixX<Scalar> F_values(n_virt, n_occ);
        for (int a = 0; a < n_virt; a++) {
            for (int i = 0; i < n_occ; i++) {
                F_values(a, i) = virtual_energies[a] - occupied_energies[i];
            }
        }
        return F_values;
    }


    /**
     *  @return All the spin-orbital energies, with the alpha spin-orbital energies appearing before the beta spin-orbital energies.
     */
    VectorX<Scalar> spinOrbitalEnergiesBlocked() const {

        const auto K = this->numberOfSpatialOrbitals();

        VectorX<Scalar> total_orbital_energies {2 * K};
        total_orbital_energies.head(K) = this->orbitalEnergies();
        total_orbital_energies.tail(K) = this->orbitalEnergies();

        return total_orbital_energies;
    }


    /**
     *  @return All the spin-orbital energies, with the alpha and beta-spinorbital energies interleaved.
     */
    VectorX<Scalar> spinOrbitalEnergiesInterleaved() const {

        const auto K = this->numberOfSpatialOrbitals();

        VectorX<Scalar> total_orbital_energies {2 * K};
        for (size_t p = 0; p < K; p++) {
            total_orbital_energies(2 * p) = this->orbitalEnergy(p);
            total_orbital_energies(2 * p + 1) = this->orbitalEnergy(p);
        }

        return total_orbital_energies;
    }


    /**
     *  @return The orbital energies belonging to the virtual orbitals.
     */
    std::vector<Scalar> virtualOrbitalEnergies() const {

        // Determine the number of occupied orbitals.
        const auto n_occ = this->orbitalSpace().numberOfOrbitals(OccupationType::k_occupied);

        std::vector<Scalar> mo_energies;  // We use a std::vector in order to be able to slice the vector later on.
        for (int i = 0; i < this->numberOfSpatialOrbitals(); i++) {
            mo_energies.push_back(this->orbitalEnergy(i));
        }

        // Add the values with indices greater than the occupied orbital indices, i.e. the virtual orbital indices, to the new vector.
        std::vector<Scalar> mo_energies_virtual;
        std::copy(mo_energies.begin() + n_occ, mo_energies.end(), std::back_inserter(mo_energies_virtual));
        return mo_energies_virtual;
    }


    /**
     *  @return The orbital energies belonging to the occupied orbitals.
     */
    std::vector<Scalar> occupiedOrbitalEnergies() const {

        // Determine the number of occupied orbitals.
        const auto n_occ = this->orbitalSpace().numberOfOrbitals(OccupationType::k_occupied);

        std::vector<Scalar> mo_energies;  // We use a std::vector in order to be able to slice the vector later on.
        for (int i = 0; i < this->numberOfSpatialOrbitals(); i++) {
            mo_energies.push_back(this->orbitalEnergy(i));
        }

        // Add the values with indices smaller than the occupied orbital indices, to the new vector.
        std::vector<Scalar> mo_energies_occupied;
        std::copy(mo_energies.begin(), mo_energies.begin() + n_occ, std::back_inserter(mo_energies_occupied));
        return mo_energies_occupied;
    }


    /**
     *  @return All the spatial orbital energies.
     */
    const VectorX<Scalar>& orbitalEnergies() const { return this->orbital_energies; }

    /**
     *  @param i            The index of the orbital.
     *
     *  @return The i-th orbital energy.
     */
    Scalar orbitalEnergy(const size_t i) const { return this->orbital_energies(i); }


    /*
     *  MARK: Error
     */

    /**
     *  @param F                The Fock operator expressed in a scalar basis.
     *  @param D                The RHF density matrix in the same scalar basis.
     *  @param S                The overlap operator of that scalar basis.
     *
     *  @return The RHF error matrix.
     */
    static SquareMatrix<Scalar> calculateError(const ScalarRSQOneElectronOperator<Scalar>& F, const Orbital1DM<Scalar>& D, const ScalarRSQOneElectronOperator<Scalar>& S) {
        return F.parameters() * D.matrix() * S.parameters() - S.parameters() * D.matrix() * F.parameters();
    }


    /*
     *  MARK: Hessians
     */

    /**
     *  @param sq_hamiltonian       The Hamiltonian expressed in an orthonormal basis.
     *  @param N_P                  The number of electron pairs.
     *  @param a                    The first virtual orbital index.
     *  @param i                    The first occupied orbital index.
     *  @param b                    The second virtual orbital index.
     *  @param j                    The second occupied orbital index.
     *
     *  @return An element of the RHF orbital Hessian.
     */
    static Scalar calculateOrbitalHessianElement(const RSQHamiltonian<Scalar>& sq_hamiltonian, const size_t N_P, const size_t a, const size_t i, const size_t b, const size_t j) {

        // Prepare some variables.
        const auto& g_op = sq_hamiltonian.twoElectron();
        const auto K = g_op.numberOfOrbitals();  // The number of spatial orbitals.
        const auto& g = g_op.parameters();

        const auto orbital_space = RHF<Scalar>::orbitalSpace(K, N_P);


        double value {0.0};

        // Inactive Fock matrix part.
        const auto F = sq_hamiltonian.calculateInactiveFockian(orbital_space).parameters();
        if (i == j) {
            value += F(a, b);
        }

        if (a == b) {
            value -= F(i, j);
        }


        // The two-electron part.
        value += 4 * g(a, i, b, j) - g(a, b, i, j) - g(a, j, b, i);

        return 4 * value;
    }


    /**
     *  @param sq_hamiltonian       The Hamiltonian expressed in an orthonormal basis.
     *  @param N_P                  The number of electron pairs.
     *
     *  @return The RHF orbital Hessian as a ImplicitRankFourTensorSlice, i.e. an object with a suitable operator() implemented.
     */
    static ImplicitRankFourTensorSlice<Scalar> calculateOrbitalHessianTensor(const RSQHamiltonian<Scalar>& sq_hamiltonian, const size_t N_P) {

        // Create an occupied-virtual orbital space.
        const auto K = sq_hamiltonian.numberOfOrbitals();
        const auto orbital_space = OrbitalSpace::Implicit({{OccupationType::k_occupied, N_P}, {OccupationType::k_virtual, K - N_P}});  // N_P occupied spatial orbitals, K-N_P virtual spatial orbitals

        auto hessian = orbital_space.template initializeRepresentableObjectFor<Scalar>(OccupationType::k_virtual, OccupationType::k_occupied, OccupationType::k_virtual, OccupationType::k_occupied);  // zero-initialize an object suitable for the representation of a virtual-occupied,virtual-occupied object (ai,bj)

        // Loop over all indices (ai,bj) to construct the orbital hessian.
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
     *  Calculate the RHF orbital Hessian (H_RI -i H_II), which can be used as a response force constant when solving the CP(R)HF equations for a purely imaginary response.
     *
     *  @param sq_hamiltonian               the Hamiltonian expressed in the canonical RHF orbital basis, resulting from a real optimization
     *  @param orbital_space                the orbital space that encapsulates the occupied-virtual separation
     *
     *  @return An RHF orbital Hessian.
     */
    template <typename Z = Scalar>
    static enable_if_t<std::is_same<Z, double>::value, ImplicitRankFourTensorSlice<complex>> calculateOrbitalHessianForImaginaryResponse(const RSQHamiltonian<double>& sq_hamiltonian, const OrbitalSpace& orbital_space) {

        using namespace GQCP::literals;

        // Prepare some variables.
        const auto& g = sq_hamiltonian.twoElectron().parameters();
        const auto F = sq_hamiltonian.calculateInactiveFockian(orbital_space).parameters();


        // Zero-initialize a virtual-occupied-virtual-occupied object.
        auto H = orbital_space.initializeRepresentableObjectFor<complex>(OccupationType::k_virtual, OccupationType::k_occupied, OccupationType::k_virtual, OccupationType::k_occupied);


        // Calculate the elements of the specific orbital Hessian (ai,bj)
        for (const auto& a : orbital_space.indices(OccupationType::k_virtual)) {
            for (const auto& i : orbital_space.indices(OccupationType::k_occupied)) {
                for (const auto& b : orbital_space.indices(OccupationType::k_virtual)) {
                    for (const auto& j : orbital_space.indices(OccupationType::k_occupied)) {

                        complex value {};  // The term in parentheses.

                        // Add the contribution from the orbital energies, i.e. the diagonal term.
                        if ((a == b) && (i == j)) {
                            value += F(a, a) - F(i, i);
                        }

                        // Add the contributions from the other terms.
                        value += g(i, b, j, a) - g(i, j, b, a);

                        H(a, i, b, j) = -4_ii * value;
                    }
                }
            }
        }

        return H;
    }


    /*
     *  MARK: Stability
     */

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
    MatrixX<Scalar> calculateSingletAStabilityMatrix(const RSQHamiltonian<Scalar>& rsq_hamiltonian) const {

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
                        singlet_A_slice(i, a, j, b) = 2.0 * g(a, i, j, b) - g(a, b, j, i);
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
    MatrixX<Scalar> calculateSingletBStabilityMatrix(const RSQHamiltonian<Scalar>& rsq_hamiltonian) const {

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
                        singlet_B_slice(i, a, j, b) = 2.0 * g(a, i, b, j) - g(a, j, b, i);
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
    MatrixX<Scalar> calculateTripletAStabilityMatrix(const RSQHamiltonian<Scalar>& rsq_hamiltonian) const {

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
    MatrixX<Scalar> calculateTripletBStabilityMatrix(const RSQHamiltonian<Scalar>& rsq_hamiltonian) const {

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


    /*
     *  MARK: Density matrices
     */

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

        // The 1-DM for RHF looks like (for K=5, N=6):
        //    2  0  0  0  0
        //    0  2  0  0  0
        //    0  0  2  0  0
        //    0  0  0  0  0
        //    0  0  0  0  0

        SquareMatrix<Scalar> D = SquareMatrix<Scalar>::Zero(K);
        D.topLeftCorner(N / 2, N / 2) = 2 * SquareMatrix<Scalar>::Identity(N / 2);

        return Orbital1DM<Scalar> {D};
    }


    /**
     *  @return The 1-DM expressed in an orthonormal spinor basis related to these optimal RHF parameters.
     */
    Orbital1DM<Scalar> calculateOrthonormalBasis1DM() const {

        const auto K = this->numberOfSpatialOrbitals();
        const auto N = 2 * this->numberOfElectronPairs();
        return RHF<Scalar>::calculateOrthonormalBasis1DM(K, N);
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
            throw std::invalid_argument("QCMethod::RHF::calculateOrthonormalBasis2DM(const size_t, const size_t): The number of given electrons cannot be odd for RHF.");
        }

        const size_t N_P = N / 2;  // The number of electron pairs.


        // Create the orbital space to determine the loops.
        const auto orbital_space = RHF<Scalar>::orbitalSpace(K, N_P);


        // Implement a KISS formula for the RHF 2-DM.
        SquareRankFourTensor<Scalar> d = SquareRankFourTensor<Scalar>::Zero(K);

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

        return Orbital2DM<Scalar>(d);
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
     *  @param C    The coefficient matrix that expresses every spatial orbital (as a column) in its underlying scalar basis.
     *  @param N    The number of electrons.
     *
     *  @return The RHF 1-DM expressed in the underlying scalar basis.
     */
    static Orbital1DM<Scalar> calculateScalarBasis1DM(const RTransformation<Scalar>& C, const size_t N) {

        const size_t K = C.numberOfOrbitals();
        const auto D_orthonormal = RHF<Scalar>::calculateOrthonormalBasis1DM(K, N);

        // Transform the 1-DM in an orthonormal basis to the underlying scalar basis.
        return D_orthonormal.transformed(C.inverse());
    }


    /**
     *  @return The RHF 1-DM in the scalar/AO basis related to these optimal RHF parameters.
     */
    Orbital1DM<Scalar> calculateScalarBasis1DM() const {

        const auto N = 2 * this->numberOfElectronPairs();
        return RHF<Scalar>::calculateScalarBasis1DM(this->expansion(), N);
    }


    /**
     *  @param C    The coefficient matrix that expresses every spatial orbital (as a column) in its underlying scalar basis.
     *  @param N    The number of electrons.
     *
     *  @return The RHF 2-DM expressed in the underlying scalar basis.
     */
    static Orbital2DM<Scalar> calculateScalarBasis2DM(const RTransformation<Scalar>& C, const size_t N) {

        const size_t K = C.numberOfOrbitals();
        const auto d_orthonormal = RHF<Scalar>::calculateOrthonormalBasis2DM(K, N);

        // Transform the 1-DM in an orthonormal basis to the underlying scalar basis.
        return d_orthonormal.transformed(C.inverse());
    }


    /**
     *  @return The RHF 1-DM in the scalar/AO basis related to these optimal RHF parameters.
     */
    Orbital2DM<Scalar> calculateScalarBasis2DM() const {

        const auto N = 2 * this->numberOfElectronPairs();
        return RHF<Scalar>::calculateScalarBasis2DM(this->expansion(), N);
    }


    /**
     *  Calculate the RHF Fock operator F = H_core + G, in which G is a contraction of the density matrix and the two-electron integrals.
     *
     *  @param D                    The RHF density matrix in a scalar basis.
     *  @param sq_hamiltonian       The Hamiltonian expressed in the same scalar basis.
     *
     *  @return The RHF Fock operator expressed in the scalar basis.
     */
    static ScalarRSQOneElectronOperator<Scalar> calculateScalarBasisFockMatrix(const Orbital1DM<Scalar>& D, const RSQHamiltonian<Scalar>& sq_hamiltonian) {

        // Get the two-electron parameters.
        const auto& g = sq_hamiltonian.twoElectron().parameters();

        // To calculate G, we must perform two double contractions:
        //      1. (mu nu|rho lambda) P(lambda rho),
        const Tensor<Scalar, 2> direct_contraction = g.template einsum<2>("ijkl,kl->ij", D.matrix());
        //      2. -0.5 (mu lambda|rho nu) P(lambda rho).
        const Tensor<Scalar, 2> exchange_contraction = -0.5 * g.template einsum<2>("ijkl,kj->il", D.matrix());

        // The previous contractions are Tensor<Scalar, 2> instances. In order to calculate the total G matrix, we will convert them back into GQCP::Matrix<Scalar>.
        auto G1 = direct_contraction.asMatrix();
        auto G2 = exchange_contraction.asMatrix();

        return ScalarRSQOneElectronOperator<Scalar> {sq_hamiltonian.core().parameters() + G1 + G2};
    }


    /*
     *  MARK: Orbital indices
     */


    /**
     *  @param N            The number of electrons.
     *
     *  @return The (spatial orbital, not spin-orbital) index of the RHF HOMO in an implicit orbital space.
     */
    static size_t homoIndex(const size_t N) {

        if (N % 2 != 0) {
            throw std::invalid_argument("QCModel::RHF::homoIndex(const size_t): Can't calculate the RHF HOMO index for an odd number of electrons N.");
        }

        return N / 2 - 1;  // We need to subtract 1 because computer indices start at 0.
    }


    /**
     *  @param N            The number of electrons.
     *
     *  @return The (spatial orbital, not spin-orbital) index of the RHF HOMO in an implicit orbital space.
     */
    size_t homoIndex() const { return RHF<Scalar>::homoIndex(this->numberOfElectrons()); }


    /**
     *  @param K            The number of spatial orbitals.
     *  @param N            The number of electrons.
     *
     *  @return The (spatial orbital, not spin-orbital) index of the RHF LUMO in an implicit orbital space.
     */
    static size_t lumoIndex(const size_t K, const size_t N) {

        if (N >= 2 * K) {
            throw std::invalid_argument("QCModel::RHF::lumoIndex(const size_t, constsize_t): There is no LUMO for the given number of electrons N and spatial orbitals K");
        }

        return RHF<Scalar>::homoIndex(N) + 1;
    }


    /**
     *  @param K            The number of spatial orbitals.
     *  @param N            The number of electrons.
     *
     *  @return The (spatial orbital, not spin-orbital) index of the RHF LUMO in an implicit orbital space.
     */
    size_t lumoIndex() const { return RHF<Scalar>::lumoIndex(this->numberOfSpatialOrbitals(), this->numberOfElectrons()); }


    /**
     *  @param K            The number of spatial orbitals.
     *  @param N_P          The number of electrons.
     *
     *  @return The implicit (i.e. with ascending and contiguous orbital indices) occupied-virtual orbital space that corresponds to these RHF model parameters.
     */
    static OrbitalSpace orbitalSpace(const size_t K, const size_t N_P) {

        return OrbitalSpace::Implicit({{OccupationType::k_occupied, N_P}, {OccupationType::k_virtual, K - N_P}});
    }


    /**
     *  @return The implicit occupied-virtual orbital space that is associated to these RHF model parameters.
     */
    OrbitalSpace orbitalSpace() const { return RHF<Scalar>::orbitalSpace(this->numberOfSpatialOrbitals(), this->numberOfElectronPairs()); }


    /*
     *  MARK: General information
     */

    /**
     *  @return The transformation that expresses the RHF MOs in terms of the underlying AOs.
     */
    const RTransformation<Scalar>& expansion() const { return this->C; }

    /**
     *  @return The number of electron pairs that these RHF model parameters describe.
     */
    size_t numberOfElectronPairs() const { return this->N_P; }

    /**
     *  @return The total number of electrons that these RHF model parameters describe.
     */
    size_t numberOfElectrons() const { return 2 * this->numberOfElectronPairs(); }

    /**
     *  @param sigma            The spin of the electrons (alpha or beta).
     *
     *  @return The number of sigma-electrons that these RHF model parameters describe.
     */
    size_t numberOfElectrons(const Spin sigma) const { return this->numberOfElectronPairs(); }

    /**
     *  @return The number of spatial orbitals that these RHF model parameters describe.
     */
    size_t numberOfSpatialOrbitals() const { return this->expansion().numberOfOrbitals(); }


    /*
     *  MARK: Response forces
     */

    /**
     *  Calculate the RHF response force for the perturbation due to a magnetic field.
     *
     *  @param L_op             The angular momentum operator expressed in the RHF orbital basis.
     *
     *  @return The RHF response force for the perturbation due to a magnetic field. Every column of the returned matrix contains the response force along the corresponding component: x, y, z.
     */
    template <typename Z = Scalar>
    enable_if_t<std::is_same<Z, double>::value, Matrix<complex, Dynamic, 3>> calculateMagneticFieldResponseForce(const VectorRSQOneElectronOperator<complex>& L_op) const {

        // Prepare some variables.
        const auto orbital_space = this->orbitalSpace();
        const auto& L = L_op.allParameters();

        // Every column of the matrix `F_kappa_B` contains the response force along the given component: x, y, z.
        const auto dim = orbital_space.numberOfExcitations(OccupationType::k_occupied, OccupationType::k_virtual);
        Matrix<complex, Dynamic, 3> F_kappa_B = Matrix<complex, Dynamic, 3>::Zero(dim, 3);
        for (size_t m = 0; m < 3; m++) {  // `m` labels a Cartesian direction.

            // Initialize a virtual-occupied object for every component.
            auto F_kappa_B_m = orbital_space.template initializeRepresentableObjectFor<complex>(OccupationType::k_virtual, OccupationType::k_occupied);

            for (const auto& a : orbital_space.indices(OccupationType::k_virtual)) {
                for (const auto& i : orbital_space.indices(OccupationType::k_occupied)) {
                    F_kappa_B_m(a, i) = -2.0 * L[m](i, a);
                }
            }
            F_kappa_B.col(m) = F_kappa_B_m.asVector();
        }

        return F_kappa_B;
    }


    /**
     *  Calculate the RHF response force for the perturbation due to a gauge origin translation of the magnetic field.
     *
     *  @param p_op             The linear momentum operator expressed in the RHF orbital basis.
     *
     *  @return The RHF response force for the perturbation due to a gauge origin translation of the magnetic field. Every column of the returned matrix contains the response force along the corresponding component: xy, xz, yx, yz, zx, zy.
     */
    template <typename Z = Scalar>
    enable_if_t<std::is_same<Z, double>::value, Matrix<complex, Dynamic, 6>> calculateGaugeOriginTranslationResponseForce(const VectorRSQOneElectronOperator<complex>& p_op) const {

        // Prepare some variables.
        const LeviCivitaTensor<double> epsilon {};
        const auto& p = p_op.allParameters();
        const auto orbital_space = this->orbitalSpace();


        // Every column of the matrix `F_kappa_G_mn` contains the response force along the given component: xy, xz, yx, yz, zx, zy.
        const auto dim = orbital_space.numberOfExcitations(OccupationType::k_occupied, OccupationType::k_virtual);
        Matrix<complex, Dynamic, 6> F_kappa_G = Matrix<complex, Dynamic, 6>::Zero(dim, 6);
        size_t column_index = 0;
        for (size_t m = 0; m < 3; m++) {      // `m` labels a Cartesian direction.
            for (size_t n = 0; n < 3; n++) {  // `n` labels another Cartesian direction.
                if (m == n) {                 // Skip the diagonal components: xx, yy, zz.
                    continue;
                }

                // Initialize a virtual-occupied object for every component.
                auto F_kappa_G_mn = orbital_space.template initializeRepresentableObjectFor<complex>(OccupationType::k_virtual, OccupationType::k_occupied);

                for (const auto& a : orbital_space.indices(GQCP::OccupationType::k_virtual)) {
                    for (const auto& i : orbital_space.indices(GQCP::OccupationType::k_occupied)) {
                        const auto f = epsilon.nonZeroIndex(m, n);

                        F_kappa_G_mn(a, i) = 2.0 * epsilon(m, n, f) * p[f](i, a);
                    }
                }
                F_kappa_G.col(column_index) = F_kappa_G_mn.asVector();
                column_index++;
            }
        }

        return F_kappa_G;
    }


    /*
     *  MARK: Response properties
     */

    /**
     *  Calculate the magnetic inducibility on the given point using the ipsocentric CSGT method.
     *
     *  @param r                        The point.
     *  @param orbital_space            The RHF occupied-virtual separation.
     *  @param x_B                      The RHF linear response with respect to the magnetic field perturbation.
     *  @param x_G                      The RHF linear response with respect to the gauge origin translation perturbation.
     *  @param j_op                     The second-quantized representation of the field-free current density operator in the RHF orbital basis.
     *
     *  @return The magnetic inducibility evaluated at the given point.
     */
    static Matrix<complex, 3, 3> calculateIpsocentricMagneticInducibility(const Vector<double, 3>& r, const OrbitalSpace& orbital_space, const Matrix<complex, Dynamic, 3>& x_B, const Matrix<complex, Dynamic, 6>& x_G, const VectorEvaluableRSQOneElectronOperator<CurrentDensityMatrixElement<complex, CartesianGTO>>& j_op) {

        using namespace GQCP::literals;

        // Prepare some variables.
        const auto& j = j_op.allParameters();
        Matrix<complex, 3, 3> J = Matrix<complex, 3, 3>::Zero();

        // Loop over both components of the magnetic inducibility.
        for (size_t u = 0; u < 3; u++) {      // `u` loops over the component of the induced current.
            for (size_t m = 0; m < 3; m++) {  // `m` loops over the component of the applied external magnetic field.

                // Initialize a more useful matrix representation of the linear response coefficients related to the magnetic field perturbation.
                auto x_m_matrix = MatrixX<complex>::FromColumnMajorVector(x_B.col(m), orbital_space.numberOfOrbitals(OccupationType::k_virtual), orbital_space.numberOfOrbitals(OccupationType::k_occupied));
                auto x_m = orbital_space.createRepresentableObjectFor<complex>(OccupationType::k_virtual, OccupationType::k_occupied, x_m_matrix);

                for (const auto& a : orbital_space.indices(OccupationType::k_virtual)) {
                    for (const auto& i : orbital_space.indices(OccupationType::k_occupied)) {

                        // Determine the left factor, i.e. the contribution from the field-free current density matrix elements.
                        complex left_value = 2.0_ii * j[u](i, a)(r);

                        // Determine the right factor, i.e. the contribution from the linear response coefficients.

                        // 1. The linear response coefficient related to the magnetic field perturbation.
                        complex right_value = -1.0_ii * x_m(a, i);

                        // 2. The linear response coefficient related to the gauge origin translation perturbation.
                        for (size_t n = 0; n < 3; n++) {

                            // 'm' cannot be equal to 'n': there is no such linear response.
                            if (n == m) {
                                continue;
                            }

                            // Determine the compound index 'mn' of the dyadic Cartesian direction in the 6-column matrix representation.
                            const auto row_major_index = 3 * m + n;
                            const auto mn = row_major_index < 4 ? row_major_index - 1 : row_major_index - 2;

                            // Initialize a more useful matrix representation of the linear response coefficients related to the gauge origin translation perturbation.
                            auto epsilon_mn_matrix = MatrixX<complex>::FromColumnMajorVector(x_G.col(mn), orbital_space.numberOfOrbitals(OccupationType::k_virtual), orbital_space.numberOfOrbitals(OccupationType::k_occupied));
                            auto epsilon_mn = orbital_space.createRepresentableObjectFor<complex>(OccupationType::k_virtual, OccupationType::k_occupied, epsilon_mn_matrix);

                            // Use the ipsocentric CSGT step, i.e. d = r - G_0 with G_0 the gauge origin (i.e. the origin).
                            right_value += -1.0_ii * epsilon_mn(a, i) * r(n);
                        }

                        J(u, m) += 2_ii * left_value * right_value;
                        assert(J(u, m).imag() < 1.0e-12);  // The magnetic inducibility should be a real-valued quantity.
                    }
                }
            }
        }

        return J;
    }


    /**
     *  Calculate the magnetic inducibility on the given grid using the ipsocentric CSGT method.
     *
     *  @param grid                     The grid.
     *  @param orbital_space            The RHF occupied-virtual separation.
     *  @param x_B                      The RHF linear response with respect to the magnetic field perturbation.
     *  @param x_G                      The RHF linear response with respect to the gauge origin translation perturbation.
     *  @param j_op                     The second-quantized representation of the field-free current density operator in the RHF orbital basis.
     *
     *  @return The magnetic inducibility evaluated on the given grid.
     */
    static Field<Matrix<complex, 3, 3>> calculateIpsocentricMagneticInducibility(const CubicGrid& grid, const OrbitalSpace& orbital_space, const Matrix<complex, Dynamic, 3>& x_B, const Matrix<complex, Dynamic, 6>& x_G, const VectorEvaluableRSQOneElectronOperator<CurrentDensityMatrixElement<complex, CartesianGTO>>& j_op) {

        using namespace GQCP::literals;

        // Evaluate the value of the magnetic inducibility for every point of the grid.
        std::vector<Matrix<complex, 3, 3>> J_field_values;
        J_field_values.reserve(grid.numberOfPoints());

        grid.forEach([&orbital_space, &x_B, &x_G, &J_field_values, &j_op](const Vector<double, 3>& r) {
            const auto J = RHF<complex>::calculateIpsocentricMagneticInducibility(r, orbital_space, x_B, x_G, j_op);
            J_field_values.push_back(J);
        });

        return Field<Matrix<complex, 3, 3>> {J_field_values};
    }
};


}  // namespace QCModel
}  // namespace GQCP
