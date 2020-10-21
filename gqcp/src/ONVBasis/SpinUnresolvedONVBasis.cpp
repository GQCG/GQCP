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

#include "ONVBasis/SpinUnresolvedONVBasis.hpp"

#include <boost/math/special_functions.hpp>
#include <boost/numeric/conversion/converter.hpp>


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param M            the number of spinors
 *  @param N            the number of electrons
 */
SpinUnresolvedONVBasis::SpinUnresolvedONVBasis(const size_t M, const size_t N) :
    BaseONVBasis(M, SpinUnresolvedONVBasis::calculateDimension(M, N)),
    ONVManipulator(N) {

    // Create a zero matrix of dimensions (M+1)x(N+1)
    this->vertex_weights = std::vector<std::vector<size_t>>(this->M + 1, std::vector<size_t>(this->N + 1, 0));

    // M=5   N=2
    // [ 0 0 0 ]
    // [ 0 0 0 ]
    // [ 0 0 0 ]
    // [ 0 0 0 ]
    // [ 0 0 0 ]
    // [ 0 0 0 ]


    // The largest (reverse lexical) string is the one that includes the first (M-N+1) vertices of the first column
    //      This is because every vertical move from (p,m) to (p+1,m+1) corresponds to "orbital p+1 is unoccupied".
    //      Therefore, the largest reverse lexical string is the one where the first (M-N) orbitals are unoccupied.
    //      This means that there should be (M-N) vertical moves from (0,0).
    // Therefore, we may only set the weights of first (M-N+1) vertices of the first column to 1.
    for (size_t p = 0; p < this->M - this->N + 1; p++) {
        this->vertex_weights[p][0] = 1;
    }

    // M=5   N=2
    // [ 1 0 0 ]
    // [ 1 0 0 ]
    // [ 1 0 0 ]
    // [ 1 0 0 ]
    // [ 0 0 0 ]
    // [ 0 0 0 ]


    // The recurrence relation for the vertex weights is as follows:
    //      Every element is the sum of the values of the element vertically above and the element left diagonally above.
    //      W(p,m) = W(p-1,m) + W(p-1,m-1)

    for (size_t m = 1; m < this->N + 1; m++) {
        for (size_t p = m; p < (this->M - this->N + m) + 1; p++) {
            this->vertex_weights[p][m] = this->vertex_weights[p - 1][m] + this->vertex_weights[p - 1][m - 1];
        }
    }

    // M=5   N=2
    // [ 1 0 0 ]
    // [ 1 1 0 ]
    // [ 1 2 1 ]
    // [ 1 3 3 ]
    // [ 0 4 6 ]
    // [ 0 0 10]
}


/*
 *  STATIC PUBLIC METHODS
 */

/**
 *  @param M        the number of orbitals
 *  @param N        the number of electrons
 *
 *  @return the dimension of this basis
 */
size_t SpinUnresolvedONVBasis::calculateDimension(const size_t M, const size_t N) {
    auto dim_double = boost::math::binomial_coefficient<double>(static_cast<unsigned>(M), static_cast<unsigned>(N));
    try {
        return boost::numeric::converter<size_t, double>::convert(dim_double);
    } catch (boost::numeric::bad_numeric_cast& e) {
        throw std::overflow_error("SpinUnresolvedONVBasis::calculateDimension(size_t, size_t): " + std::string(e.what()));
    }
}


/*
 *  PUBLIC OVERRIDDEN METHODS
 */

/**
 *  @param representation      a representation of an ONV
 *
 *  @return the address (i.e. the ordering number) of the given ONV
 */
size_t SpinUnresolvedONVBasis::addressOf(const size_t representation) const {

    // An implementation of the formula in Helgaker, starting the addressing count from zero
    size_t copy = representation;
    size_t address = 0;
    size_t electron_count = 0;            // counts the number of electrons in the spin string up to orbital p
    while (copy != 0) {                   // we will remove the least significant bit each loop, we are finished when no bits are left
        size_t p = __builtin_ctzl(copy);  // p is the orbital index counter (starting from 1)
        electron_count++;                 // each bit is an electron hence we add it up to the electron count
        address += this->vertexWeight(p, electron_count);
        copy ^= copy & -copy;  // flip the least significant bit
    }
    return address;
}


/**
 *  @param onv       the ONV
 *
 *  @return the number of ONVs (with a larger address) the ONV would couple with given a one electron operator
 */
size_t SpinUnresolvedONVBasis::countOneElectronCouplings(const SpinUnresolvedONV& onv) const {

    size_t V = this->M - N;  // number of virtual orbitals
    size_t coupling_count = 0;

    for (size_t e1 = 0; e1 < this->N; e1++) {
        size_t p = onv.occupationIndexOf(e1);
        coupling_count += (V + e1 - p);  // number of virtuals with an index larger than p
    }

    return coupling_count;
}


/**
 *  @return the number of non-zero couplings of a one electron coupling scheme in this ONV basis
 */
size_t SpinUnresolvedONVBasis::countTotalOneElectronCouplings() const {
    return (this->M - N) * N * (dim);
}


/**
 *  @return the number of non-zero couplings of a two electron coupling scheme in this ONV basis
 */
size_t SpinUnresolvedONVBasis::countTotalTwoElectronCouplings() const {

    size_t two_electron_permutation = 0;  // all distributions for two electrons over the virtual orbitals
    if (this->M - N >= 2) {
        two_electron_permutation = calculateDimension(this->M - N, 2) * N * (N - 1) * (dim) / 2;
    }

    return two_electron_permutation + countTotalOneElectronCouplings();
}


/**
 *  @param onv       the ONV
 *
 *  @return the number of ONVs (with a larger address) the ONV would couple with given a two electron operator
 */
size_t SpinUnresolvedONVBasis::countTwoElectronCouplings(const SpinUnresolvedONV& onv) const {

    size_t V = this->M - N;  // number of virtual orbitals
    size_t coupling_count = 0;

    for (size_t e1 = 0; e1 < this->N; e1++) {

        size_t p = onv.occupationIndexOf(e1);
        coupling_count += (V + e1 - p);  //  one electron part

        for (size_t e2 = e1 + 1; e2 < this->N; e2++) {

            size_t q = onv.occupationIndexOf(e2);
            size_t coupling_count2 = (V + e2 - q);
            coupling_count += (V - coupling_count2) * coupling_count2;

            if (coupling_count2 > 1) {
                coupling_count += calculateDimension(coupling_count2, 2);
            }
        }
    }

    return coupling_count;
}


/**
 *  Evaluate the operator in a dense matrix
 *
 *  @param one_op               the one-electron operator in an orthonormal orbital basis to be evaluated in this ONV basis
 *  @param diagonal_values      bool to indicate if diagonal values will be calculated
 *
 *  @return the operator's evaluation in a dense matrix with the dimensions of this ONV basis
 */
SquareMatrix<double> SpinUnresolvedONVBasis::evaluateOperatorDense(const ScalarGSQOneElectronOperator<double>& one_op, const bool diagonal_values) const {

    auto K = one_op.numberOfOrbitals();
    if (K != this->M) {
        throw std::invalid_argument("SpinUnresolvedONVBasis::evaluateOperatorDense(ScalarGSQOneElectronOperator<double>, bool): Basis functions of this ONV basis and the operator are incompatible.");
    }

    MatrixRepresentationEvaluationContainer<SquareMatrix<double>> evaluation_iterator(this->dim);
    this->evaluate<SquareMatrix<double>>(one_op, evaluation_iterator, diagonal_values);
    return evaluation_iterator.evaluation();
}


/**
 *  Evaluate the operator in a dense matrix
 *
 *  @param two_op               the two-electron operator in an orthonormal orbital basis to be evaluated in this ONV basis
 *  @param diagonal_values      bool to indicate if diagonal values will be calculated
 *
 *  @return the operator's evaluation in a dense matrix with the dimensions of this ONV basis
 */
SquareMatrix<double> SpinUnresolvedONVBasis::evaluateOperatorDense(const ScalarGSQTwoElectronOperator<double>& two_op, const bool diagonal_values) const {

    const auto K = two_op.numberOfOrbitals();
    if (K != this->M) {
        throw std::invalid_argument("SpinUnresolvedONVBasis::evaluateOperatorDense(ScalarGSQTwoElectronOperator<double>, bool): Basis functions of this ONV basis and the operator are incompatible.");
    }

    MatrixRepresentationEvaluationContainer<SquareMatrix<double>> evaluation_iterator(this->dim);
    this->evaluate<SquareMatrix<double>>(two_op, evaluation_iterator, diagonal_values);
    return evaluation_iterator.evaluation();
}


/**
 *  Evaluate the Hamiltonian in a dense matrix
 *
 *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal basis
 *  @param diagonal_values              bool to indicate if diagonal values will be calculated
 *
 *  @return the Hamiltonian's evaluation in a dense matrix with the dimensions of this ONV basis
 */
SquareMatrix<double> SpinUnresolvedONVBasis::evaluateOperatorDense(const RSQHamiltonian<double>& sq_hamiltonian, const bool diagonal_values) const {

    auto K = sq_hamiltonian.numberOfOrbitals();
    if (K != this->M) {
        throw std::invalid_argument("SpinUnresolvedONVBasis::evaluateOperatorDense(RSQHamiltonian<double>, bool): Basis functions of this ONV basis and the operator are incompatible.");
    }

    MatrixRepresentationEvaluationContainer<SquareMatrix<double>> evaluation_iterator(this->dim);
    this->evaluate<SquareMatrix<double>>(sq_hamiltonian.core(), sq_hamiltonian.twoElectron(), evaluation_iterator, diagonal_values);
    return evaluation_iterator.evaluation();
}


/**
 *  Evaluate the diagonal of the operator
 *
 *  @param one_op               the one-electron operator in an orthonormal orbital basis to be evaluated in this ONV basis
 *
 *  @return the operator's diagonal evaluation in a vector with the dimension of this ONV basis
 */
VectorX<double> SpinUnresolvedONVBasis::evaluateOperatorDiagonal(const ScalarGSQOneElectronOperator<double>& one_op) const {

    const auto K = one_op.numberOfOrbitals();

    if (K != this->M) {
        throw std::invalid_argument("SpinUnresolvedONVBasis::evaluateOperatorDiagonal(ScalarGSQOneElectronOperator<double>): Basis functions of this ONV basis and the operator are incompatible.");
    }

    const auto& one_op_par = one_op.parameters();

    size_t N = this->numberOfElectrons();
    size_t dim = this->dimension();

    VectorX<double> diagonal = VectorX<double>::Zero(dim);

    SpinUnresolvedONV onv = this->constructONVFromAddress(0);  // onv with address 0
    for (size_t I = 0; I < dim; I++) {                         // I loops over all addresses in this ONV basis

        if (I > 0) {
            this->transformONVToNextPermutation(onv);
        }

        for (size_t e1 = 0; e1 < N; e1++) {  // A1 (annihilation 1)
            size_t p = onv.occupationIndexOf(e1);
            diagonal(I) += one_op_par(p, p);
        }
    }

    return diagonal;
};


/**
 *  Evaluate the diagonal of the operator
 *
 *  @param two_op               the two-electron operator in an orthonormal orbital basis to be evaluated in this ONV basis
 *
 *  @return the operator's diagonal evaluation in a vector with the dimension of this ONV basis
 */
VectorX<double> SpinUnresolvedONVBasis::evaluateOperatorDiagonal(const ScalarGSQTwoElectronOperator<double>& two_op) const {

    auto K = two_op.numberOfOrbitals();
    if (K != this->M) {
        throw std::invalid_argument("SpinUnresolvedONVBasis::evaluateOperatorDiagonal(ScalarGSQTwoElectronOperator<double>): Basis functions of this ONV basis and the operator are incompatible.");
    }

    size_t N = this->numberOfElectrons();
    size_t dim = this->dimension();

    VectorX<double> diagonal = VectorX<double>::Zero(dim);

    const auto k_par = two_op.effectiveOneElectronPartition().parameters();
    const auto& two_op_par = two_op.parameters();

    SpinUnresolvedONV onv = this->constructONVFromAddress(0);  // onv with address 0
    for (size_t I = 0; I < dim; I++) {                         // I loops over all addresses in this ONV basis

        if (I > 0) {
            this->transformONVToNextPermutation(onv);
        }

        for (size_t e1 = 0; e1 < N; e1++) {  // A1 (annihilation 1)
            size_t p = onv.occupationIndexOf(e1);
            diagonal(I) += k_par(p, p);

            for (size_t q = 0; q < K; q++) {  // q loops over SOs
                if (onv.isOccupied(q)) {
                    diagonal(I) += 0.5 * two_op_par(p, p, q, q);
                } else {
                    diagonal(I) += 0.5 * two_op_par(p, q, q, p);
                }
            }
        }
    }

    return diagonal;
};


/**
 *  Evaluate the diagonal of the Hamiltonian
 *
 *  @param sq_hamiltonian           the Hamiltonian expressed in an orthonormal basis
 *
 *  @return the Hamiltonian's diagonal evaluation in a vector with the dimension of this ONV basis
 */
VectorX<double> SpinUnresolvedONVBasis::evaluateOperatorDiagonal(const RSQHamiltonian<double>& sq_hamiltonian) const {

    auto K = sq_hamiltonian.numberOfOrbitals();
    if (K != this->M) {
        throw std::invalid_argument("SpinUnresolvedONVBasis::evaluateOperatorDiagonal(RSQHamiltonian<double>): Basis functions of this ONV basis and the operator are incompatible.");
    }

    return this->evaluateOperatorDiagonal(sq_hamiltonian.core()) + this->evaluateOperatorDiagonal(sq_hamiltonian.twoElectron());
};


/**
 *  Evaluate the operator in a sparse matrix
 *
 *  @param one_op               the one-electron operator in an orthonormal orbital basis to be evaluated in this ONV basis
 *  @param diagonal_values      bool to indicate if diagonal values will be calculated
 *
 *  @return the operator's evaluation in a sparse matrix with the dimensions of this ONV basis
 */
Eigen::SparseMatrix<double> SpinUnresolvedONVBasis::evaluateOperatorSparse(const ScalarGSQOneElectronOperator<double>& one_op, const bool diagonal_values) const {

    auto K = one_op.numberOfOrbitals();
    if (K != this->M) {
        throw std::invalid_argument("SpinUnresolvedONVBasis::evaluateOperatorSparse(ScalarGSQOneElectronOperator<double>, bool): Basis functions of this ONV basis and the operator are incompatible.");
    }

    MatrixRepresentationEvaluationContainer<Eigen::SparseMatrix<double>> evaluation_iterator(this->dim);

    size_t memory = this->countTotalOneElectronCouplings();
    if (diagonal_values) {
        memory += this->dim;
    }

    evaluation_iterator.reserve(memory);
    this->evaluate<Eigen::SparseMatrix<double>>(one_op, evaluation_iterator, diagonal_values);
    evaluation_iterator.addToMatrix();
    return evaluation_iterator.evaluation();
}


/**
 *  Evaluate the operator in a sparse matrix
 *
 *  @param two_op               the two-electron operator in an orthonormal orbital basis to be evaluated in this ONV basis
 *  @param diagonal_values      bool to indicate if diagonal values will be calculated
 *
 *  @return the operator's evaluation in a sparse matrix with the dimensions of this ONV basis
 */
Eigen::SparseMatrix<double> SpinUnresolvedONVBasis::evaluateOperatorSparse(const ScalarGSQTwoElectronOperator<double>& two_op, const bool diagonal_values) const {

    auto K = two_op.numberOfOrbitals();
    if (K != this->M) {
        throw std::invalid_argument("SpinUnresolvedONVBasis::evaluateOperatorSparse(ScalarGSQTwoElectronOperator<double>, bool): Basis functions of this ONV basis and the operator are incompatible.");
    }

    MatrixRepresentationEvaluationContainer<Eigen::SparseMatrix<double>> evaluation_iterator(this->dim);

    size_t memory = this->countTotalTwoElectronCouplings();
    if (diagonal_values) {
        memory += this->dim;
    }

    evaluation_iterator.reserve(memory);
    this->evaluate<Eigen::SparseMatrix<double>>(two_op, evaluation_iterator, diagonal_values);
    evaluation_iterator.addToMatrix();
    return evaluation_iterator.evaluation();
}


/**
 *  Evaluate the Hamiltonian in a sparse matrix
 *
 *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal basis
 *  @param diagonal_values              bool to indicate if diagonal values will be calculated
 *
 *  @return the Hamiltonian's evaluation in a sparse matrix with the dimensions of this ONV basis
 */
Eigen::SparseMatrix<double> SpinUnresolvedONVBasis::evaluateOperatorSparse(const RSQHamiltonian<double>& sq_hamiltonian, const bool diagonal_values) const {

    auto K = sq_hamiltonian.numberOfOrbitals();
    if (K != this->M) {
        throw std::invalid_argument("SpinUnresolvedONVBasis::evaluateOperatorSparse(RSQHamiltonian<double>, bool): Basis functions of this ONV basis and the operator are incompatible.");
    }

    MatrixRepresentationEvaluationContainer<Eigen::SparseMatrix<double>> evaluation_iterator(this->dim);

    size_t memory = this->countTotalTwoElectronCouplings();
    if (diagonal_values) {
        memory += this->dim;
    }

    evaluation_iterator.reserve(memory);
    this->evaluate<Eigen::SparseMatrix<double>>(sq_hamiltonian.core(), sq_hamiltonian.twoElectron(), evaluation_iterator, diagonal_values);
    evaluation_iterator.addToMatrix();
    return evaluation_iterator.evaluation();
}


/**
 *  Calculate unsigned representation for a given address
 *
 *  @param address                 the address of the representation is calculated
 *
 *  @return unsigned representation of the address
 */
size_t SpinUnresolvedONVBasis::representationOf(size_t address) const {
    size_t representation = 0;
    if (this->N != 0) {
        representation = 0;
        size_t m = this->N;  // counts the number of electrons in the spin string up to orbital p

        for (size_t p = this->M; p > 0; p--) {  // p is an orbital index
            size_t weight = vertexWeight(p - 1, m);

            if (weight <= address) {  // the algorithm can move diagonally, so we found an occupied orbital
                address -= weight;
                representation |= ((1) << (p - 1));  // set the (p-1)th bit: see (https://stackoverflow.com/a/47990)

                m--;  // since we found an occupied orbital, we have one electron less
                if (m == 0) {
                    break;
                }
            }
        }
    }
    return representation;
}


/**
 *  @param representation       a representation of an ONV
 *
 *  @return the next bitstring permutation
 *
 *      Examples:
 *          011 -> 101
 *          101 -> 110
 */
size_t SpinUnresolvedONVBasis::nextPermutationOf(const size_t representation) const {

    // t gets this->representation's least significant 0 bits set to 1
    unsigned long t = representation | (representation - 1UL);

    // Next set to 1 the most significant bit to change,
    // set to 0 the least significant ones, and add the necessary 1 bits.
    return (t + 1UL) | (((~t & (t + 1UL)) - 1UL) >> (__builtin_ctzl(representation) + 1UL));
}


/*
 *  PUBLIC METHODS
 */

/**
 *  Calculates sigma(pq) + sigma(qp)'s: all one-electron couplings for each annihilation-creation pair in this ONV basis
 *  and stores them in sparse matrices for each pair combination
 *
 *  @return a vector of sparse matrices containing the one-electron couplings for the spin-unresolved basis
 *      Ordered as: sigma(00), sigma(01) + sigma(10), sigma(02)+ sigma(20), ...
 */
std::vector<Eigen::SparseMatrix<double>> SpinUnresolvedONVBasis::calculateOneElectronCouplings() const {

    size_t K = this->numberOfOrbitals();
    size_t N = this->numberOfElectrons();
    size_t dim = this->dimension();

    std::vector<std::vector<Eigen::Triplet<double>>> sparse_entries(K * (K + 1) / 2);
    std::vector<Eigen::SparseMatrix<double>> sparse_matrices(K * (K + 1) / 2, Eigen::SparseMatrix<double>(dim, dim));

    if (N == 0) {
        return sparse_matrices;
    }

    // Reserve appropriate number of entries
    size_t reservation_size = SpinUnresolvedONVBasis::calculateDimension(K - 1, N - 1);
    for (size_t p = 0; p < K; p++) {
        sparse_entries[p * (K + K + 1 - p) / 2].reserve(reservation_size);
        for (size_t q = p + 1; q < K; q++) {
            sparse_entries[p * (K + K + 1 - p) / 2 + q - p].reserve(2 * reservation_size);
        }
    }

    SpinUnresolvedONV onv = this->constructONVFromAddress(0);  // onv with address 0
    for (size_t I = 0; I < dim; I++) {                         // I loops over all the addresses of the onv
        for (size_t e1 = 0; e1 < N; e1++) {                    // e1 (electron 1) loops over the (number of) electrons
            size_t p = onv.occupationIndexOf(e1);              // retrieve the index of a given electron
            // remove the weight from the initial address I, because we annihilate
            size_t address = I - this->vertexWeight(p, e1 + 1);
            // The e2 iteration counts the number of encountered electrons for the creation operator
            // We only consider greater addresses than the initial one (because of symmetry)
            // Hence we only count electron after the annihilated electron (e1)
            sparse_entries[p * (K + K + 1 - p) / 2].emplace_back(I, I, 1);
            size_t e2 = e1 + 1;
            size_t q = p + 1;
            int sign_e2 = 1;
            // perform a shift
            this->shiftUntilNextUnoccupiedOrbital<1>(onv, address, q, e2, sign_e2);
            size_t dex = 0;
            while (q < K) {
                size_t J = address + this->vertexWeight(q, e2);

                sparse_entries[p * (K + K + 1 - p) / 2 + q - p].emplace_back(I, J, sign_e2);
                sparse_entries[p * (K + K + 1 - p) / 2 + q - p].emplace_back(J, I, sign_e2);

                q++;  // go to the next orbital
                dex++;
                // perform a shift
                this->shiftUntilNextUnoccupiedOrbital<1>(onv, address, q, e2, sign_e2);
            }  //  (creation)
        }      // e1 loop (annihilation)

        // Prevent last permutation
        if (I < dim - 1) {
            this->transformONVToNextPermutation(onv);
        }
    }

    for (size_t k = 0; k < K * (K + 1) / 2; k++) {
        sparse_matrices[k].setFromTriplets(sparse_entries[k].begin(), sparse_entries[k].end());
    }

    return sparse_matrices;
}


/**
 *  Evaluate a one electron operator in a matrix vector product
 *
 *  @param one_op                       the one electron operator expressed in an orthonormal basis
 *  @param x                            the vector upon which the evaluation acts 
 *  @param diagonal                     the diagonal evaluated in this ONV basis
 *
 *  @return the one electron operator's matrix vector product in a vector with the dimensions of this ONV basis
 */
VectorX<double> SpinUnresolvedONVBasis::evaluateOperatorMatrixVectorProduct(const ScalarGSQOneElectronOperator<double>& one_op, const VectorX<double>& x, const VectorX<double>& diagonal) const {

    auto K = one_op.numberOfOrbitals();
    if (K != this->M) {
        throw std::invalid_argument("SpinUnresolvedONVBasis::evaluateOperatorMatrixVectorProduct(ScalarGSQOneElectronOperator<double>, VectorX<double>, VectorX<double>): Basis functions of this ONV basis and the operator are incompatible.");
    }

    MatrixRepresentationEvaluationContainer<VectorX<double>> evaluation_iterator(x, diagonal);
    this->evaluate<VectorX<double>>(one_op, evaluation_iterator, false);
    return evaluation_iterator.evaluation();
}


/**
 *  Evaluate a two electron operator in a matrix vector product
 *
 *  @param two_op                       the two electron operator expressed in an orthonormal basis
 *  @param x                            the vector upon which the evaluation acts 
 *  @param diagonal                     the diagonal evaluated in this ONV basis
 *
 *  @return the two electron operator's matrix vector product in a vector with the dimensions of this ONV basis
 */
VectorX<double> SpinUnresolvedONVBasis::evaluateOperatorMatrixVectorProduct(const ScalarGSQTwoElectronOperator<double>& two_op, const VectorX<double>& x, const VectorX<double>& diagonal) const {

    auto K = two_op.numberOfOrbitals();
    if (K != this->M) {
        throw std::invalid_argument("SpinUnresolvedONVBasis::evaluateOperatorMatrixVectorProduct(ScalarGSQTwoElectronOperator<double>, VectorX<double>, VectorX<double>): Basis functions of this ONV basis and the operator are incompatible.");
    }

    MatrixRepresentationEvaluationContainer<VectorX<double>> evaluation_iterator(x, diagonal);
    this->evaluate<VectorX<double>>(two_op, evaluation_iterator, false);
    return evaluation_iterator.evaluation();
}


/**
 *  Evaluate the Hamiltonian in a matrix vector product
 *
 *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal basis
 *  @param x                            the vector upon which the evaluation acts 
 *  @param diagonal                     the diagonal evaluated in this ONV basis
 *
 *  @return the Hamiltonian's matrix vector product in a vector with the dimensions of this ONV basis
 */
VectorX<double> SpinUnresolvedONVBasis::evaluateOperatorMatrixVectorProduct(const RSQHamiltonian<double>& sq_hamiltonian, const VectorX<double>& x, const VectorX<double>& diagonal) const {

    auto K = sq_hamiltonian.numberOfOrbitals();
    if (K != this->M) {
        throw std::invalid_argument("SpinUnresolvedONVBasis::evaluateOperatorMatrixVectorProduct(RSQHamiltonian<double>, VectorX<double>, VectorX<double>): Basis functions of this ONV basis and the operator are incompatible.");
    }

    MatrixRepresentationEvaluationContainer<VectorX<double>> evaluation_iterator(x, diagonal);
    this->evaluate<VectorX<double>>(sq_hamiltonian.core(), sq_hamiltonian.twoElectron(), evaluation_iterator, false);
    return evaluation_iterator.evaluation();
}


/**
 *  Iterate over all ONVs in this ONV basis and apply the given callback function.
 * 
 *  @param callback             the function to be applied in every iteration. Its arguments are a spin-unresolved ONV and its corresponding addresses
 */
void SpinUnresolvedONVBasis::forEach(const std::function<void(const SpinUnresolvedONV&, const size_t)>& callback) const {

    SpinUnresolvedONV onv = this->constructONVFromAddress(0);

    for (size_t I = 0; I < this->dim; I++) {  // I loops over addresses of all ONVs

        callback(onv, I);

        if (I < this->dim - 1) {  // prevent the last permutation from occurring
            this->transformONVToNextPermutation(onv);
        }
    }  // address (I) loop
}


}  // namespace GQCP
