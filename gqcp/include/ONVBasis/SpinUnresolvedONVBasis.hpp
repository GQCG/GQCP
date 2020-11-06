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


#include "Mathematical/Representation/MatrixRepresentationEvaluationContainer.hpp"
#include "ONVBasis/ONVPath.hpp"
#include "ONVBasis/SpinUnresolvedONV.hpp"
#include "Operator/SecondQuantized/GSQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/USQOneElectronOperatorComponent.hpp"

#include <functional>


namespace GQCP {


/**
 *  The full spin-unresolved ONV basis for a number of spinors/spin-orbitals and number of electrons.
 */
class SpinUnresolvedONVBasis {
private:
    // The number of spinors/spin-orbitals.
    size_t M;

    // The number of electrons, i.e. the number of occupied spinors/spin-orbitals.
    size_t N;

    // The vertex weights corresponding to the addressing scheme for a full spin-unresolved ONV basis. This addressing scheme is taken from Helgaker, Jørgensen, Olsen (2000).
    std::vector<std::vector<size_t>> vertex_weights;

public:
    // The ONV that is naturally related to a full spin-unresolved ONV basis. See also `ONVPath`.
    using ONV = SpinUnresolvedONV;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  @param M            The number of spinors/spin-orbitals.
     *  @param N            The number of electrons, i.e. the number of occupied spinors/spin-orbitals.
     */
    SpinUnresolvedONVBasis(const size_t M, const size_t N);

    /**
     *  The default constructor.
     */
    SpinUnresolvedONVBasis() = default;


    /*
     *  MARK: Basic information
     */

    /**
     *  Calculate the dimension of the full spin-unresolved ONV basis with a given number of spinors/spin-orbitals and number of electrons.
     * 
     *  @param M            The number of spinors/spin-orbitals.
     *  @param N            The number of electrons, i.e. the number of occupied spinors/spin-orbitals.
     *
     *  @return The dimension of a spin-unresolved ONV basis.
     */
    static size_t calculateDimension(const size_t M, const size_t N);

    /**
     *  @return The dimension of this ONV basis.
     */
    size_t dimension() const { return SpinUnresolvedONVBasis::calculateDimension(this->numberOfOrbitals(), this->numberOfElectrons()); }

    /**
     *  @return The number of electrons, i.e. the number of occupied spinors/spin-orbitals.
     */
    size_t numberOfElectrons() const { return this->N; }

    /**
     *  @return The number of spinors/spin-orbitals.
     */
    size_t numberOfOrbitals() const { return this->M; }


    /**
     *  MARK: Addressing scheme, address calculations and ONV manipulations
     */

    /**
     *  Access the arc weight of an arc in the addressing scheme of this ONV basis. The addressing scheme is taken from Helgaker, Jørgensen, Olsen (2000).
     * 
     *  @param p            The orbital index.
     *  @param n            The electron index.
     *
     *  @return The arc weight of the arc starting at the given vertex (p, n).
     */
    size_t arcWeight(const size_t p, const size_t n) const;

    /**
     *  @param p            The orbital index.
     *  @param n            The electron index.
     * 
     *  @return The vertex weight related to the given indices (p,n).
     */
    size_t vertexWeight(const size_t p, const size_t n) const { return this->vertex_weights[p][n]; }

    /**
     *  @return All the vertex weights for this ONV basis, stored as a vector of vectors. The outer axis represents the orbital indices, the inner axis represents the electron indices.
     */
    const std::vector<std::vector<size_t>>& vertexWeights() const { return this->vertex_weights; }

    /**
     *  Calculate the address (i.e. the ordering number) of an unsigned representation of a spin-unresolved ONV.
     * 
     *  @param representation      The unsigned representation of a spin-unresolved ONV.
     *
     *  @return The address corresponding to the unsigned representation of a spin-unresolved ONV.
     */
    size_t addressOf(const size_t representation) const;

    /**
     *  Calculate the address (i.e. the ordering number) of a spin-unresolved ONV.
     * 
     *  @param onv          The spin-unresolved ONV.
     *
     *  @return The address (i.e. the ordering number) of the given spin-unresolved ONV.
     */
    size_t addressOf(const SpinUnresolvedONV& onv) const { return this->addressOf(onv.unsignedRepresentation()); }

    /**
     *  Calculate the next allowed unsigned representation of a spin-unresolved ONV in this ONV basis.
     * 
     *  @param representation       A representation of a spin-unresolved ONV.
     *
     *  @return The next allowed unsigned representation of a spin-unresolved ONV in this ONV basis.
     *
     *  @example
     *          011 -> 101
     *          101 -> 110
     */
    size_t nextPermutationOf(const size_t representation) const;

    /**
     *  Calculate the unsigned representation of a spin-unresolved ONV that corresponds to the address/ordering number in this ONV basis.
     *
     *  @param address                 The address/ordering number of a spin-unresolved ONV in this ONV basis.
     *
     *  @return The unsigned representation of a spin-unresolved ONV that corresponds to the address/ordering number in this ONV basis.
     */
    size_t representationOf(const size_t address) const;

    /**
     *  Create the ONV that corresponds to the given address in this ONV basis.
     * 
     *  @param address                 The address/ordering number of a spin-unresolved ONV in this ONV basis.
     *
     *  @return The ONV that corresponds to the given address in this ONV basis.
     */
    SpinUnresolvedONV constructONVFromAddress(const size_t address) const;

    /**
     *  Modify a `SpinResolvedONV` to the next allowed ONV in this ONV basis.
     * 
     *  @param onv      A spin-resolved ONV.
     */
    void transformONVToNextPermutation(SpinUnresolvedONV& onv) const;

    /**
     *  Modify a `SpinResolvedONV` to the one with the given address in this ONV basis.
     *
     *  @param onv          A spin-resolved ONV.
     *  @param address      The target address in this ONV basis.
     */
    void transformONVCorrespondingToAddress(SpinUnresolvedONV& onv, const size_t address) const;


    /*
     *  MARK: Couplings
     */

    /**
     *  Calculate the number of ONVs (with a larger address) that a given spin-unresolved ONV would couple with in this ONV basis, through a one-electron operator.
     * 
     *  @param onv          The spin-unresolved ONV.
     *
     *  @return The number of ONVs (with a larger address) that a given spin-unresolved ONV would couple with in this ONV basis, through a one-electron operator.
     */
    size_t countOneElectronCouplings(const SpinUnresolvedONV& onv) const;

    /**
     *  @return The total number of non-zero and non-diagonal couplings of a one-electron operator in this ONV basis.
     */
    size_t countTotalOneElectronCouplings() const;

    /**
     *  @return The total number of non-zero and non-diagonal couplings of a two-electron operator in this ONV basis.
     */
    size_t countTotalTwoElectronCouplings() const;

    /**
     *  Calculate the number of ONVs (with a larger address) that a given spin-unresolved ONV would couple with in this ONV basis, through a two-electron operator.
     * 
     *  @param onv          The spin-unresolved ONV.
     *
     *  @return The number of ONVs (with a larger address) that a given spin-unresolved ONV would couple with in this ONV basis, through a two-electron operator.
     */
    size_t countTwoElectronCouplings(const SpinUnresolvedONV& onv) const;

    /**
     *  Calculate all one-electron coupling elements for this spin-unresolved ONV basis. These are all the intermediate matrices sigma(pq)_{IJ}, as defined by Helgaker, Jørgensen, Olsen (2000).
     *
     *  @return A vector of sparse matrices containing the one-electron coupling elements for this spin-unresolved ONV basis. The elements of this vector are ordered through the one-electron excitation (pq) in ascending order: sigma(00), sigma(01) + sigma(10), sigma(02)+ sigma(20), ...
     */
    std::vector<Eigen::SparseMatrix<double>> calculateOneElectronCouplings() const;


    /**
     *  MARK: Iterating
     */

    /**
     *  Iterate over all ONVs in this ONV basis and apply the given callback function.
     * 
     *  @param callback            The function to be applied in every iteration. Its supplied arguments are a spin-unresolved ONV and its corresponding addresses.
     */
    void forEach(const std::function<void(const SpinUnresolvedONV&, const size_t)>& callback) const;


    /*
     *  MARK: Dense operator evaluations
     */

    /**
     *  Calculate the dense matrix representation of a generalized one-electron operator in this ONV basis.
     *
     *  @param f                A generalized one-electron operator expressed in an orthonormal orbital basis.
     *
     *  @return A dense matrix represention of the one-electron operator.
     */
    SquareMatrix<double> evaluateOperatorDense(const ScalarGSQOneElectronOperator<double>& f) const;

    /**
     *  Calculate the dense matrix representation of a component of an unrestricted one-electron operator in this ONV basis.
     *
     *  @param f                A a component of an unrestricted one-electron operator expressed in an orthonormal orbital basis.
     *
     *  @return A dense matrix represention of the one-electron operator.
     */
    SquareMatrix<double> evaluateOperatorDense(const ScalarUSQOneElectronOperatorComponent<double>& f) const {

        // We may convert an unrestricted component into the generalized representation.
        const auto f_generalized = ScalarGSQOneElectronOperator<double>::FromUnrestrictedComponent(f);
        return this->evaluateOperatorDense(f_generalized);
    }


    /*
     *  MARK: Diagonal restricted operator evaluations
     */


    /*
     *  MARK: Sparse restricted operator evaluations
     */

    /*
     *  MARK: Restricted matrix-vector product evaluations.
     */


    /*
     *  MARK: Operator evaluations - general implementations - containers
     */

    /**
     *  Calculate the matrix representation of a generalized one-electron operator in this ONV basis and emplace it in the given container.
     * 
     *  @tparam Matrix                      The type of matrix used to store the evaluations.
     *
     *  @param f_op                         A generalized one-electron operator expressed in an orthonormal spinor basis.
     *  @param container                    A specialized container for emplacing evaluations/matrix elements.
     */
    template <typename Matrix>
    void evaluate(const ScalarGSQOneElectronOperator<double>& f_op, MatrixRepresentationEvaluationContainer<Matrix>& container) const {

        const auto& f = f_op.parameters();
        const auto dim = this->dimension();

        SpinUnresolvedONV onv = this->constructONVFromAddress(0);  // start with ONV with address 0

        for (; !container.isFinished(); container.increment()) {  // loops over all possible ONVs
            for (size_t e1 = 0; e1 < N; e1++) {                   // loop over electrons that can be annihilated

                // Create an ONVPath for each new ONV.
                ONVPath<SpinUnresolvedONVBasis> onv_path {*this, onv};

                size_t q = onv.occupationIndexOf(e1);  // retrieve orbital index of the electron that will be annihilated

                // The diagonal values are a result of annihilation-creation on the same orbital index and are thus the same as the initial ONV.
                container.addRowwise(container.index, f(q, q));

                // For the non-diagonal values, we will create all possible matrix elements of the Hamiltonian in the routine below.
                onv_path.annihilate(q, e1);

                // Stop the loop if 1) the path is finished, meaning that orbital index p is at M (the total number of orbitals) and 2) if the orbital index is out of bounds after left translation of a vertical arc.
                while (!onv_path.isFinished() && onv_path.isOrbitalIndexValid()) {

                    // Find the next unoccupied orbital, i.e. the next vertical arc in the path.
                    onv_path.leftTranslateDiagonalArcUntilVerticalArc();

                    // Calculate the address of the path if we would close it right now.
                    const size_t address = onv_path.addressAfterCreation();

                    const double value = onv_path.sign() * f(onv_path.orbitalIndex(), q);

                    // Add the one-electron integral as matrix elements of a Hermitian matrix.
                    container.addColumnwise(address, value);
                    container.addRowwise(address, value);

                    // Move orbital index such that other unoccupied orbitals can be found within the loop.
                    onv_path.leftTranslateVerticalArc();
                }
            }
            // Prevent last ONV since there is no possibility for an electron to be annihilated anymore.
            if (container.index < dim - 1) {
                this->transformONVToNextPermutation(onv);
            }
        }
    }


    // /**
    //  *  Evaluate the operator in a dense matrix
    //  *
    //  *  @param two_op               the two-electron operator in an orthonormal orbital basis to be evaluated in the spin-unresolved ONV basis
    //  *  @param diagonal_values      bool to indicate if diagonal values will be calculated
    //  *
    //  *  @return the operator's evaluation in a dense matrix with the dimensions of the spin-unresolved ONV basis
    //  */
    // SquareMatrix<double> evaluateOperatorDense(const ScalarGSQTwoElectronOperator<double>& two_op, const bool diagonal_values) const;

    // /**
    //  *  Evaluate the Hamiltonian in a dense matrix
    //  *
    //  *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal basis
    //  *  @param diagonal_values              bool to indicate if diagonal values will be calculated
    //  *
    //  *  @return the Hamiltonian's evaluation in a dense matrix with the dimensions of the spin-unresolved ONV basis
    //  */
    // SquareMatrix<double> evaluateOperatorDense(const GSQHamiltonian<double>& sq_hamiltonian, const bool diagonal_values) const;

    // /**
    //  *  Evaluate the diagonal of the operator
    //  *
    //  *  @param f_op               the one-electron operator in an orthonormal orbital basis to be evaluated in the spin-unresolved ONV basis
    //  *
    //  *  @return the operator's diagonal evaluation in a vector with the dimension of the spin-unresolved ONV basis
    //  */
    // VectorX<double> evaluateOperatorDiagonal(const ScalarGSQOneElectronOperator<double>& f_op) const;

    // /**
    //  *  Evaluate the diagonal of the operator
    //  *
    //  *  @param two_op               the two-electron operator in an orthonormal orbital basis to be evaluated in the spin-unresolved ONV basis
    //  *
    //  *  @return the operator's diagonal evaluation in a vector with the dimension of the spin-unresolved ONV basis
    //  */
    // VectorX<double> evaluateOperatorDiagonal(const ScalarGSQTwoElectronOperator<double>& two_op) const;

    // /**
    //  *  Evaluate the diagonal of the Hamiltonian
    //  *
    //  *  @param sq_hamiltonian           the Hamiltonian expressed in an orthonormal basis
    //  *
    //  *  @return the Hamiltonian's diagonal evaluation in a vector with the dimension of the spin-unresolved ONV basis
    //  */
    // VectorX<double> evaluateOperatorDiagonal(const GSQHamiltonian<double>& sq_hamiltonian) const;

    // /**
    //  *  Evaluate the operator in a sparse matrix
    //  *
    //  *  @param f_op               the one-electron operator in an orthonormal orbital basis to be evaluated in the spin-unresolved ONV basis
    //  *  @param diagonal_values      bool to indicate if diagonal values will be calculated
    //  *
    //  *  @return the operator's evaluation in a sparse matrix with the dimensions of the spin-unresolved ONV basis
    //  */
    // Eigen::SparseMatrix<double> evaluateOperatorSparse(const ScalarGSQOneElectronOperator<double>& f_op, const bool diagonal_values) const;

    // /**
    //  *  Evaluate the operator in a sparse matrix
    //  *
    //  *  @param two_op               the two-electron operator in an orthonormal orbital basis to be evaluated in the spin-unresolved ONV basis
    //  *  @param diagonal_values      bool to indicate if diagonal values will be calculated
    //  *
    //  *  @return the operator's evaluation in a sparse matrix with the dimensions of the spin-unresolved ONV basis
    //  */
    // Eigen::SparseMatrix<double> evaluateOperatorSparse(const ScalarGSQTwoElectronOperator<double>& two_op, const bool diagonal_values) const;

    // /**
    //  *  Evaluate the Hamiltonian in a sparse matrix
    //  *
    //  *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal basis
    //  *  @param diagonal_values              bool to indicate if diagonal values will be calculated
    //  *
    //  *  @return the Hamiltonian's evaluation in a sparse matrix with the dimensions of the spin-unresolved ONV basis
    //  */
    // Eigen::SparseMatrix<double> evaluateOperatorSparse(const RSQHamiltonian<double>& sq_hamiltonian, const bool diagonal_values) const;


    // PUBLIC METHODS


    /**
     *  Evaluate the operator in a given evaluation iterator in the spin-unresolved ONV basis
     *
     *  @tparam Matrix                       the type of matrix used to store the evaluations
     *
     *  @param f_op                        the one-electron operator in an orthonormal orbital basis to be evaluated in the spin-unresolved ONV basis
     *  @param evaluation_iterator           evaluation iterator to which the evaluations are added
     *  @param diagonal_values               bool to indicate if diagonal values will be calculated
     */
    // template <typename _Matrix>
    // void evaluate(const ScalarGSQOneElectronOperator<double>& f_op, MatrixRepresentationEvaluationContainer<_Matrix>& evaluation_iterator, const bool diagonal_values) const {

    //     const auto& f = f_op.parameters();

    //     const size_t K = this->numberOfOrbitals();
    //     const size_t N = this->numberOfElectrons();
    //     const size_t dim = this->dimension();

    //     SpinUnresolvedONV onv = this->constructONVFromAddress(0);  // onv with address 0

    //     for (; !evaluation_iterator.isFinished(); evaluation_iterator.increment()) {  // I loops over all the addresses of the onv
    //         for (size_t e1 = 0; e1 < N; e1++) {                                       // e1 (electron 1) loops over the (number of) electrons

    //             size_t p = onv.occupationIndexOf(e1);  // retrieve the index of a given electron
    //             // remove the weight from the initial address I, because we annihilate
    //             size_t address = evaluation_iterator.index - this->vertexWeight(p, e1 + 1);

    //             if (diagonal_values) {
    //                 evaluation_iterator.addRowwise(evaluation_iterator.index, f(p, p));
    //             }

    //             // The e2 iteration counts the number of encountered electrons for the creation operator
    //             // We only consider greater addresses than the initial one (because of symmetry)
    //             // Hence we only count electron after the annihilated electron (e1)
    //             size_t e2 = e1 + 1;
    //             size_t q = p + 1;

    //             int sign_e2 = 1;
    //             // perform a shift
    //             this->shiftUntilNextUnoccupiedOrbital<1>(onv, address, q, e2, sign_e2);

    //             while (q < K) {
    //                 size_t J = address + this->vertexWeight(q, e2);
    //                 double value = sign_e2 * f(p, q);
    //                 evaluation_iterator.addColumnwise(J, value);
    //                 evaluation_iterator.addRowwise(J, value);

    //                 q++;  // go to the next orbital

    //                 // perform a shift
    //                 this->shiftUntilNextUnoccupiedOrbital<1>(onv, address, q, e2, sign_e2);
    //             }  //  (creation)
    //         }      // e1 loop (annihilation)
    //         // Prevent last permutation
    //         if (evaluation_iterator.index < dim - 1) {
    //             this->transformONVToNextPermutation(onv);
    //         }
    //     }
    // }


    /**
     *  Evaluate a one-electron operator in this spin-unresolved ONV basis.
     * 
     *  @tparam Representation              The matrix representation that is used for storing the result. Essentially, any type that can be used in MatrixRepresentationEvaluationContainer<Representation>.
     * 
     *  @param f_op                       A one-electron operator in an orthonormal orbital basis.
     *  @param should_calculate_diagonal    If diagonal values should be calculated.
     * 
     *  @return The matrix representation of the given one-electron operator.
     */
    // template <typename Representation>
    // Representation evaluate_old(const ScalarGSQOneElectronOperator<double>& f_op, const bool diagonal_values) const {

    //     const auto& f = f_op.parameters();

    //     MatrixRepresentationEvaluationContainer<Representation> evaluation_iterator {this->dimension()};

    //     const size_t K = this->numberOfOrbitals();
    //     const size_t N = this->numberOfElectrons();
    //     const size_t dim = this->dimension();

    //     SpinUnresolvedONV onv = this->constructONVFromAddress(0);  // onv with address 0

    //     for (; !evaluation_iterator.isFinished(); evaluation_iterator.increment()) {  // I loops over all the addresses of the onv
    //         for (size_t e1 = 0; e1 < N; e1++) {                                       // e1 (electron 1) loops over the (number of) electrons

    //             size_t p = onv.occupationIndexOf(e1);  // retrieve the index of a given electron
    //             // remove the weight from the initial address I, because we annihilate
    //             size_t address = evaluation_iterator.index - this->vertexWeight(p, e1 + 1);

    //             if (diagonal_values) {
    //                 evaluation_iterator.addRowwise(evaluation_iterator.index, f(p, p));
    //             }

    //             // The e2 iteration counts the number of encountered electrons for the creation operator
    //             // We only consider greater addresses than the initial one (because of symmetry)
    //             // Hence we only count electron after the annihilated electron (e1)
    //             size_t e2 = e1 + 1;
    //             size_t q = p + 1;

    //             int sign_e2 = 1;
    //             // perform a shift
    //             this->shiftUntilNextUnoccupiedOrbital<1>(onv, address, q, e2, sign_e2);

    //             while (q < K) {
    //                 size_t J = address + this->vertexWeight(q, e2);
    //                 double value = sign_e2 * f(p, q);
    //                 evaluation_iterator.addColumnwise(J, value);
    //                 evaluation_iterator.addRowwise(J, value);

    //                 q++;  // go to the next orbital

    //                 // perform a shift
    //                 this->shiftUntilNextUnoccupiedOrbital<1>(onv, address, q, e2, sign_e2);
    //             }  //  (creation)
    //         }      // e1 loop (annihilation)
    //         // Prevent last permutation
    //         if (evaluation_iterator.index < dim - 1) {
    //             this->transformONVToNextPermutation(onv);
    //         }
    //     }
    //     return evaluation_iterator.evaluation();
    // }


    /**
     *  Evaluate the operator in a given evaluation iterator in the spin-unresolved ONV basis
     *
     *  @tparam Matrix                       the type of matrix used to store the evaluations
     *
     *  @param two_op                        the two-electron operator in an orthonormal orbital basis to be evaluated in the spin-unresolved ONV basis
     *  @param evaluation_iterator           evaluation iterator to which the evaluations are added
     *  @param diagonal_values               bool to indicate if diagonal values will be calculated
     */
    // template <typename _Matrix>
    // void evaluate(const ScalarGSQTwoElectronOperator<double>& two_op, MatrixRepresentationEvaluationContainer<_Matrix>& evaluation_iterator, const bool diagonal_values) const {

    //     // Calling this combined method for both the one- and two-electron operator does not affect the performance, hence we avoid writing more code by plugging a zero operator in the combined method
    //     evaluate(ScalarGSQOneElectronOperator<double> {this->M}, two_op, evaluation_iterator, diagonal_values);
    // }


    /**
     *  Evaluate the operators in a given evaluation iterator in the spin-unresolved ONV basis
     *
     *  @tparam Matrix                       the type of matrix used to store the evaluations
     *
     *  @param f_op                        the one-electron operator in an orthonormal orbital basis to be evaluated in the spin-unresolved ONV basis
     *  @param two_op                        the two-electron operator in an orthonormal orbital basis to be evaluated in the spin-unresolved ONV basis
     *  @param evaluation_iterator           evaluation iterator to which the evaluations are added
     *  @param diagonal_values               bool to indicate if diagonal values will be calculated
     */
    // template <typename _Matrix>
    // void evaluate(const ScalarGSQOneElectronOperator<double>& f_op, const ScalarGSQTwoElectronOperator<double>& two_op, MatrixRepresentationEvaluationContainer<_Matrix>& evaluation_iterator, const bool diagonal_values) const {

    //     const auto& two_op_par = two_op.parameters();

    //     const size_t K = this->numberOfOrbitals();
    //     const size_t N = this->numberOfElectrons();
    //     const size_t dim = this->dimension();

    //     ScalarGSQOneElectronOperator<double> k = two_op.effectiveOneElectronPartition() + f_op;
    //     const auto& k_par = k.parameters();

    //     SpinUnresolvedONV onv = this->constructONVFromAddress(0);  // onv with address 0

    //     for (; !evaluation_iterator.isFinished(); evaluation_iterator.increment()) {  // I loops over all addresses in the spin-unresolved ONV basis
    //         if (evaluation_iterator.index > 0) {
    //             this->transformONVToNextPermutation(onv);
    //         }
    //         int sign1 = -1;                      // start with -1 because we flip at the start of the annihilation (so we start at 1, followed by:  -1, 1, ...)
    //         for (size_t e1 = 0; e1 < N; e1++) {  // A1 (annihilation 1)

    //             sign1 *= -1;
    //             size_t p = onv.occupationIndexOf(e1);  // retrieve the index of a given electron
    //             size_t address = evaluation_iterator.index - this->vertexWeight(p, e1 + 1);

    //             // Strictly diagonal values
    //             if (diagonal_values) {
    //                 evaluation_iterator.addRowwise(evaluation_iterator.index, k_par(p, p));
    //                 for (size_t q = 0; q < K; q++) {  // q loops over SOs
    //                     if (onv.isOccupied(q)) {
    //                         evaluation_iterator.addRowwise(evaluation_iterator.index, 0.5 * two_op_par(p, p, q, q));
    //                     } else {
    //                         evaluation_iterator.addRowwise(evaluation_iterator.index, 0.5 * two_op_par(p, q, q, p));
    //                     }
    //                 }
    //             }

    //             size_t address1 = address;
    //             size_t e2 = e1;
    //             size_t q = p;

    //             /**
    //              *  A1 > C1 (annihlation 1 > creation 1)
    //              */
    //             int sign2 = sign1;
    //             q--;
    //             e2--;
    //             this->shiftUntilPreviousUnoccupiedOrbital<1>(onv, address1, q, e2, sign2);
    //             while (q != -1) {

    //                 size_t address2 = address1 + this->vertexWeight(q, e2 + 2);

    //                 /**
    //                  *  C2 > A2
    //                  */
    //                 int sign3 = sign1;
    //                 for (size_t e3 = e1 + 1; e3 < N; e3++) {
    //                     sign3 *= -1;  // initial sign3 = sign of the annihilation, with one extra electron(from crea) = *-1
    //                     size_t r = onv.occupationIndexOf(e3);
    //                     size_t address3 = address2 - this->vertexWeight(r, e3 + 1);

    //                     size_t e4 = e3 + 1;
    //                     size_t s = r + 1;

    //                     int sign4 = sign3;
    //                     this->shiftUntilNextUnoccupiedOrbital<1>(onv, address3, s, e4, sign4);

    //                     while (s < K) {
    //                         size_t J = address3 + this->vertexWeight(s, e4);
    //                         int signev = sign1 * sign2 * sign3 * sign4;
    //                         double value = signev * 0.5 * (two_op_par(p, q, r, s) + two_op_par(r, s, p, q) - two_op_par(p, s, r, q) - two_op_par(r, q, p, s));


    //                         evaluation_iterator.addColumnwise(J, value);
    //                         evaluation_iterator.addRowwise(J, value);

    //                         s++;
    //                         this->shiftUntilNextUnoccupiedOrbital<1>(onv, address3, s, e4, sign4);
    //                     }
    //                 }
    //                 q--;
    //                 this->shiftUntilPreviousUnoccupiedOrbital<1>(onv, address1, q, e2, sign2);
    //             }

    //             e2 = e1 + 1;
    //             q = p + 1;
    //             sign2 = sign1;
    //             this->shiftUntilNextUnoccupiedOrbital<1>(onv, address, q, e2, sign2);

    //             /**
    //              *  A1 < C1
    //              */
    //             while (q < K) {

    //                 address1 = address + this->vertexWeight(q, e2);

    //                 /**
    //                  *  A2 > C1
    //                  */
    //                 int sign3 = sign2;
    //                 for (size_t e3 = e2; e3 < N; e3++) {
    //                     sign3 *= -1;  // -1 cause we created electron (creation) sign of A is now the that of C *-1
    //                     size_t r = onv.occupationIndexOf(e3);
    //                     size_t address3 = address1 - this->vertexWeight(r, e3 + 1);

    //                     size_t e4 = e3 + 1;
    //                     size_t s = r + 1;
    //                     int sign4 = sign3;
    //                     this->shiftUntilNextUnoccupiedOrbital<1>(onv, address3, s, e4, sign4);

    //                     while (s < K) {
    //                         size_t J = address3 + this->vertexWeight(s, e4);
    //                         int signev = sign1 * sign2 * sign3 * sign4;

    //                         double value = signev * 0.5 * (two_op_par(p, q, r, s) + two_op_par(r, s, p, q) - two_op_par(r, q, p, s) - two_op_par(p, s, r, q));

    //                         evaluation_iterator.addColumnwise(J, value);
    //                         evaluation_iterator.addRowwise(J, value);

    //                         s++;  // go to the next orbital
    //                         this->shiftUntilNextUnoccupiedOrbital<1>(onv, address3, s, e4, sign4);

    //                     }  // (creation)
    //                 }

    //                 size_t r = q;
    //                 sign3 = sign2;
    //                 size_t address1c = address1;

    //                 /**
    //                  *  A2 < C1, (A2 > A1)
    //                  */
    //                 for (size_t e3 = e2 - 1; e3 > e1; e3--) {
    //                     sign3 *= -1;
    //                     size_t e4 = e2;
    //                     address1c += this->vertexWeight(r, e3) - this->vertexWeight(r, e3 + 1);
    //                     r = onv.occupationIndexOf(e3);
    //                     size_t address2 = address1c - this->vertexWeight(r, e3);
    //                     int sign4 = sign2;
    //                     size_t s = q + 1;
    //                     this->shiftUntilNextUnoccupiedOrbital<1>(onv, address2, s, e4, sign4);
    //                     while (s < K) {

    //                         size_t J = address2 + this->vertexWeight(s, e4);

    //                         int signev = sign1 * sign2 * sign3 * sign4;
    //                         double value = signev * 0.5 * (two_op_par(p, q, r, s) + two_op_par(r, s, p, q) - two_op_par(r, q, p, s) - two_op_par(p, s, r, q));

    //                         evaluation_iterator.addColumnwise(J, value);
    //                         evaluation_iterator.addRowwise(J, value);

    //                         s++;
    //                         this->shiftUntilNextUnoccupiedOrbital<1>(onv, address2, s, e4, sign4);
    //                     }
    //                 }

    //                 /**
    //                  *  A2 = C1
    //                  */
    //                 int signev = sign2 * sign1;

    //                 double value_I = k_par(p, q);  // cover the one electron calculations

    //                 for (size_t s = 0; s < K; s++) {
    //                     if (!onv.isOccupied(s)) {
    //                         value_I += 0.5 * (two_op_par(p, s, s, q));
    //                     } else {
    //                         value_I += 0.5 * (two_op_par(s, s, p, q) - two_op_par(s, q, p, s) + two_op_par(p, q, s, s));
    //                     }
    //                 }

    //                 value_I *= signev;

    //                 q++;

    //                 evaluation_iterator.addColumnwise(address1, value_I);
    //                 evaluation_iterator.addRowwise(address1, value_I);

    //                 this->shiftUntilNextUnoccupiedOrbital<1>(onv, address, q, e2, sign2);
    //             }
    //         }
    //     }
    // }


    /**
     *  Evaluate the matrix-vector product of a one-electron operator
     *
     *  @param f_op                       the one-electron operator expressed in an orthonormal basis
     *  @param x                            the vector of the matrix-vector product
     *  @param diagonal                     the diagonal of the matrix representation of the operator inside the spin-unresolved ONV basis
     *
     *  @return a vector that is equal to the matrix-vector product of the one-electron operator's matrix representation and the given vector
     */
    // VectorX<double> evaluateOperatorMatrixVectorProduct(const ScalarGSQOneElectronOperator<double>& f_op, const VectorX<double>& x, const VectorX<double>& diagonal) const;

    /**
     *  Evaluate a two electron operator in a matrix vector product
     *
     *  @param two_op                       the two electron operator expressed in an orthonormal basis
     *  @param x                            the vector upon which the evaluation acts 
     *  @param diagonal                     the diagonal evaluated in the spin-unresolved ONV basis
     *
     *  @return a vector that is equal to the matrix-vector product of the two-electron operator's matrix representation and the given vector
     */
    // VectorX<double> evaluateOperatorMatrixVectorProduct(const ScalarGSQTwoElectronOperator<double>& two_op, const VectorX<double>& x, const VectorX<double>& diagonal) const;

    /**
     *  Evaluate the Hamiltonian in a matrix vector product
     *
     *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal basis
     *  @param x                            the vector upon which the evaluation acts 
     *  @param diagonal                     the diagonal evaluated in the spin-unresolved ONV basis
     *
     *  @return a vector that is equal to the matrix-vector product of the Hamiltonian's matrix representation and the given vector
     */
    // VectorX<double> evaluateOperatorMatrixVectorProduct(const RSQHamiltonian<double>& sq_hamiltonian, const VectorX<double>& x, const VectorX<double>& diagonal) const;


    /*
     *  MARK: Legacy code
     */

    /**
     *  Find the next unoccupied orbital in a given spin-unresolved ONV,
     *  update the electron count, orbital index,
     *  and update the address by calculating a shift
     *  resulting from a difference between the initial vertex weights for the encountered occupied orbitals
     *  and the corrected vertex weights accounting for previously annihilated electrons
     *
     *  @tparam T        the number of previously annihilated electrons
     *
     *  @param address   the address which is updated
     *  @param onv       the spin-unresolved ONV for which we search the next unnocupied orbital
     *  @param q         the orbital index
     *  @param e         the electron count
     */
    template <int T>
    void shiftUntilNextUnoccupiedOrbital(const SpinUnresolvedONV& onv, size_t& address, size_t& q, size_t& e) const {

        // Test whether the current orbital index is occupied
        while (e < this->N && q == onv.occupationIndexOf(e)) {

            // Take the difference of vertex weights for the encountered electron weights to that of a vertex weight path with "a" fewer electrons
            // +1 is added to the electron index, because of how the addressing scheme is arrayed.
            address += this->vertexWeight(q, e + 1 - T) - this->vertexWeight(q, e + 1);

            // move to the next electron and orbital
            e++;
            q++;
        }
    }


    /**
     *  Find the next unoccupied orbital in a given spin-unresolved ONV,
     *  update the electron count, orbital index, sign,
     *  and update the address by calculating a shift
     *  resulting from a difference between the initial vertex weights for the encountered occupied orbitals
     *  and the corrected vertex weights accounting for previously annihilated electrons
     *
     *  @tparam T        the number of previously annihilated electrons
     *
     *  @param address   the address which is updated
     *  @param onv       the spin-unresolved ONV for which we search the next unnocupied orbital
     *  @param q         the orbital index
     *  @param e         the electron count
     *  @param sign      the sign which is flipped for each iteration
     */
    template <int T>
    void shiftUntilNextUnoccupiedOrbital(const SpinUnresolvedONV& onv, size_t& address, size_t& q, size_t& e, int& sign) const {

        // Test whether the current orbital index is occupied
        while (e < this->N && q == onv.occupationIndexOf(e)) {

            // Take the difference of vertex weights for the encountered electron weights to that of a vertex weight path with "a" fewer electrons
            // +1 is added to the electron index, because of how the addressing scheme is arrayed.
            address += this->vertexWeight(q, e + 1 - T) - this->vertexWeight(q, e + 1);

            // move to the next electron and orbital
            e++;
            q++;
            sign *= -1;
        }
    }


    /**
     *  Find the previous unoccupied orbital in a given spin-unresolved ONV,
     *  update the electron count, orbital index, sign,
     *  and update the address by calculating a shift
     *  resulting from a difference between the initial vertex weights for the encountered occupied orbitals
     *  and the corrected vertex weights accounting for newly created electrons
     *
     *  @tparam T        the number of newly created electrons
     *
     *  @param address   the address which is updated
     *  @param onv       the spin-unresolved ONV for which we search the next unoccupied orbital
     *  @param q         the orbital index
     *  @param e         the electron count
     *  @param sign      the sign which is flipped for each iteration
     */
    template <int T>
    void shiftUntilPreviousUnoccupiedOrbital(const SpinUnresolvedONV& onv, size_t& address, size_t& q, size_t& e, int& sign) const {

        // Test whether the current orbital index is occupied
        while (e != -1 && q == onv.occupationIndexOf(e)) {

            int shift = static_cast<int>(this->vertexWeight(q, e + 1 + T)) - static_cast<int>(this->vertexWeight(q, e + 1));
            address += shift;

            e--;
            q--;
            sign *= -1;
        }
    }
};


}  // namespace GQCP
