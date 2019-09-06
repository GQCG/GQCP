// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2019  the GQCG developers
// 
// GQCG-gqcp is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// GQCG-gqcp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-gqcp.  If not, see <http://www.gnu.org/licenses/>.
// 
#include "FockSpace/FockSpace.hpp"

#include <boost/math/special_functions.hpp>
#include <boost/numeric/conversion/converter.hpp>


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param K        the number of orbitals
 *  @param N        the number of electrons
 */
FockSpace::FockSpace(size_t K, size_t N) :
        BaseFockSpace (K, FockSpace::calculateDimension(K, N)),
        FockPermutator (N)
{
    // Create a zero matrix of dimensions (K+1)x(N+1)
    this->vertex_weights = Matrixu(this->K + 1, Vectoru(this->N + 1, 0));

    // K=5   N=2
    // [ 0 0 0 ]
    // [ 0 0 0 ]
    // [ 0 0 0 ]
    // [ 0 0 0 ]
    // [ 0 0 0 ]
    // [ 0 0 0 ]


    // The largest (reverse lexical) string is the one that includes the first (K-N+1) vertices of the first column
    //      This is because every vertical move from (p,m) to (p+1,m+1) corresponds to "orbital p+1 is unoccupied".
    //      Therefore, the largest reverse lexical string is the one where the first (K-N) orbitals are unoccupied.
    //      This means that there should be (K-N) vertical moves from (0,0).
    // Therefore, we may only set the weights of first (K-N+1) vertices of the first column to 1.
    for (size_t p = 0; p < this->K - this->N + 1; p++) {
        this->vertex_weights[p][0] = 1;
    }

    // K=5   N=2
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
        for (size_t p = m; p < (this->K - this->N + m) + 1; p++) {
            this->vertex_weights[p][m] = this->vertex_weights[p - 1][m] + this->vertex_weights[p - 1][m - 1];
        }
    }

    // K=5   N=2
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
 *  @param K        the number of orbitals
 *  @param N        the number of electrons
 *
 *  @return the dimension of the Fock space
 */
size_t FockSpace::calculateDimension(size_t K, size_t N) {
    auto dim_double = boost::math::binomial_coefficient<double>(static_cast<unsigned>(K), static_cast<unsigned>(N));
    try {
        return boost::numeric::converter<size_t, double>::convert(dim_double);
    } catch (boost::numeric::bad_numeric_cast &e) {
        throw std::overflow_error("FockSpace::calculateDimension(size_t, size_t): "+ std::string(e.what()));
    }
}



/*
 *  PUBLIC METHODS
 */

/**
 *  @param representation       a representation of an ONV
 *
 *  @return the next bitstring permutation
 *
 *      Examples:
 *          011 -> 101
 *          101 -> 110
 */
size_t FockSpace::ulongNextPermutation(size_t representation) const {

    // t gets this->representation's least significant 0 bits set to 1
    unsigned long t = representation | (representation - 1UL);

    // Next set to 1 the most significant bit to change,
    // set to 0 the least significant ones, and add the necessary 1 bits.
    return (t + 1UL) | (((~t & (t+1UL)) - 1UL) >> (__builtin_ctzl(representation) + 1UL));
}


/**
 *  @param representation      a representation of an ONV
 *
 *  @return the address (i.e. the ordering number) of the given ONV
 */
size_t FockSpace::getAddress(size_t unsigned_onv) const {
    // An implementation of the formula in Helgaker, starting the addressing count from zero
    size_t address = 0;
    size_t electron_count = 0;  // counts the number of electrons in the spin string up to orbital p
    while(unsigned_onv != 0) {  // we will remove the least significant bit each loop, we are finished when no bits are left
        size_t p = __builtin_ctzl(unsigned_onv);  // p is the orbital index counter (starting from 1)
        electron_count++;  // each bit is an electron hence we add it up to the electron count
        address += this->get_vertex_weights(p , electron_count);
        unsigned_onv ^= unsigned_onv & -unsigned_onv;  // flip the least significant bit
    }
    return address;

}


/**
 *  Calculate unsigned representation for a given address
 *
 *  @param address                 the address of the representation is calculated
 *
 *  @return unsigned representation of the address
 */
size_t FockSpace::calculateRepresentation(size_t address) const {
    size_t representation = 0;
    if (this->N != 0) {
        representation = 0;
        size_t m = this->N;  // counts the number of electrons in the spin string up to orbital p

        for (size_t p = this->K; p > 0; p--) {  // p is an orbital index
            size_t weight = get_vertex_weights(p-1, m);

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
 *  @param onv       the ONV
 *
 *  @return the amount of ONVs (with a larger address) this ONV would couple with given a one electron operator
 */
size_t FockSpace::countOneElectronCouplings(const ONV& onv) const {
    size_t V = K-N;  // amount of virtual orbitals
    size_t coupling_count = 0;

    for (size_t e1 = 0; e1 < this->N; e1++) {
        size_t p = onv.get_occupation_index(e1);
        coupling_count += (V + e1 - p);  // amount of virtuals with an index larger than p
    }

    return coupling_count;
}


/**
 *  @param onv       the ONV
 *
 *  @return the amount of ONVs (with a larger address) this ONV would couple with given a two electron operator
 */
size_t FockSpace::countTwoElectronCouplings(const ONV& onv) const {

    size_t V = K-N; // amount of virtual orbitals
    size_t coupling_count = 0;

    for (size_t e1 = 0; e1 < this->N; e1++){

        size_t p = onv.get_occupation_index(e1);
        coupling_count += (V + e1 - p);  //  one electron part

        for (size_t e2 = e1+1; e2 < this->N; e2++){

            size_t q = onv.get_occupation_index(e2);
            size_t coupling_count2 = (V + e2 - q);
            coupling_count += (V-coupling_count2)*coupling_count2;

            if(coupling_count2 > 1 ){
                coupling_count += calculateDimension(coupling_count2, 2);
            }
        }
    }

    return coupling_count;
}


/**
 *  @return the amount non-zero couplings of a one electron coupling scheme in the Fock space
 */
size_t FockSpace::countTotalOneElectronCouplings() const {
    return (K-N)*N*(dim);
}


/**
 *  @return the amount non-zero couplings of a two electron coupling scheme in the Fock space
 */
size_t FockSpace::countTotalTwoElectronCouplings() const {

    size_t two_electron_permutation = 0; // all distributions for two electrons over the virtual orbitals
    if (K-N >= 2) {
        two_electron_permutation = calculateDimension(K-N, 2)*N*(N-1)*(dim)/2;
    }

    return two_electron_permutation + countTotalOneElectronCouplings();
}


/**
 *  Evaluate the operator in a dense matrix
 *
 *  @param one_op               the one-electron operator in an orthonormal orbital basis to be evaluated in the Fock space
 *  @param diagonal_values      bool to indicate if diagonal values will be calculated
 *
 *  @return the operator's evaluation in a dense matrix with the dimensions of the Fock space
 */
SquareMatrix<double> FockSpace::evaluateOperatorDense(const ScalarSQOneElectronOperator<double>& one_op, bool diagonal_values) const {

    auto K = one_op.dimension();
    if (K != this->K) {
        throw std::invalid_argument("FockSpace::evaluateOperatorDense(ScalarSQOneElectronOperator<double>, bool): Basis functions of the Fock space and the operator are incompatible.");
    }

    EvaluationMatrix<SquareMatrix<double>> container (this->dim);
    this->EvaluateOperator<SquareMatrix<double>>(one_op, container, diagonal_values);
    return container.get_matrix();
}


/**
 *  Evaluate the operator in a sparse matrix
 *
 *  @param one_op               the one-electron operator in an orthonormal orbital basis to be evaluated in the Fock space
 *  @param diagonal_values      bool to indicate if diagonal values will be calculated
 *
 *  @return the operator's evaluation in a sparse matrix with the dimensions of the Fock space
 */
Eigen::SparseMatrix<double> FockSpace::evaluateOperatorSparse(const ScalarSQOneElectronOperator<double>& one_op, bool diagonal_values) const {

    auto K = one_op.dimension();
    if (K != this->K) {
        throw std::invalid_argument("FockSpace::evaluateOperatorSparse(ScalarSQOneElectronOperator<double>, bool): Basis functions of the Fock space and the operator are incompatible.");
    }

    EvaluationMatrix<Eigen::SparseMatrix<double>> container (this->dim);

    size_t memory = this->countTotalOneElectronCouplings();
    if (diagonal_values) {
        memory += this->dim;
    }

    container.reserve(memory);
    this->EvaluateOperator<Eigen::SparseMatrix<double>>(one_op, container, diagonal_values);
    container.addToMatrix();
    return container.get_matrix();
}


/**
 *  Evaluate the operator in a dense matrix
 *
 *  @param two_op               the two-electron operator in an orthonormal orbital basis to be evaluated in the Fock space
 *  @param diagonal_values      bool to indicate if diagonal values will be calculated
 *
 *  @return the operator's evaluation in a dense matrix with the dimensions of the Fock space
 */
SquareMatrix<double> FockSpace::evaluateOperatorDense(const ScalarSQTwoElectronOperator<double>& two_op, bool diagonal_values) const {


    const auto K = two_op.dimension();
    if (K != this->K) {
        throw std::invalid_argument("FockSpace::evaluateOperatorDense(ScalarSQTwoElectronOperator<double>, bool): Basis functions of the Fock space and the operator are incompatible.");
    }

    EvaluationMatrix<SquareMatrix<double>> container (this->dim);
    this->EvaluateOperator<SquareMatrix<double>>(two_op, container, diagonal_values);
    return container.get_matrix();
}


/**
 *  Evaluate the operator in a sparse matrix
 *
 *  @param two_op               the two-electron operator in an orthonormal orbital basis to be evaluated in the Fock space
 *  @param diagonal_values      bool to indicate if diagonal values will be calculated
 *
 *  @return the operator's evaluation in a sparse matrix with the dimensions of the Fock space
 */
Eigen::SparseMatrix<double> FockSpace::evaluateOperatorSparse(const ScalarSQTwoElectronOperator<double>& two_op, bool diagonal_values) const {

    auto K = two_op.dimension();
    if (K != this->K) {
        throw std::invalid_argument("FockSpace::evaluateOperatorSparse(ScalarSQTwoElectronOperator<double>, bool): Basis functions of the Fock space and the operator are incompatible.");
    }

    EvaluationMatrix<Eigen::SparseMatrix<double>> container (this->dim);

    size_t memory = this->countTotalTwoElectronCouplings();
    if (diagonal_values) {
        memory += this->dim;
    }

    container.reserve(memory);
    this->EvaluateOperator<Eigen::SparseMatrix<double>>(two_op, container, diagonal_values);
    container.addToMatrix();
    return container.get_matrix();
}


/**
 *  Evaluate the Hamiltonian in a dense matrix
 *
 *  @param ham_par              the Hamiltonian parameters in an orthonormal orbital basis to be evaluated in the Fock space
 *  @param diagonal_values      bool to indicate if diagonal values will be calculated
 *
 *  @return the Hamiltonian's evaluation in a dense matrix with the dimensions of the Fock space
 */
SquareMatrix<double> FockSpace::evaluateOperatorDense(const SQHamiltonian<double>& ham_par, bool diagonal_values) const {

    auto K = ham_par.dimension();
    if (K != this->K) {
        throw std::invalid_argument("FockSpace::evaluateOperatorDense(SQHamiltonian<double>, bool): Basis functions of the Fock space and the operator are incompatible.");
    }

    EvaluationMatrix<SquareMatrix<double>> container (this->dim);
    this->EvaluateOperator<SquareMatrix<double>>(ham_par.core(), ham_par.twoElectron(), container, diagonal_values);
    return container.get_matrix();
}


/**
 *  Evaluate the Hamiltonian in a sparse matrix
 *
 *  @param ham_par              the Hamiltonian parameters in an orthonormal orbital basis to be evaluated in the Fock space
 *  @param diagonal_values      bool to indicate if diagonal values will be calculated
 *
 *  @return the Hamiltonian's evaluation in a sparse matrix with the dimensions of the Fock space
 */
Eigen::SparseMatrix<double> FockSpace::evaluateOperatorSparse(const SQHamiltonian<double>& ham_par, bool diagonal_values) const {

    auto K = ham_par.dimension();
    if (K != this->K) {
        throw std::invalid_argument("FockSpace::evaluateOperatorSparse(SQHamiltonian<double>, bool): Basis functions of the Fock space and the operator are incompatible.");
    }

    EvaluationMatrix<Eigen::SparseMatrix<double>> container (this->dim);

    size_t memory = this->countTotalTwoElectronCouplings();
    if (diagonal_values) {
        memory += this->dim;
    }

    container.reserve(memory);
    this->EvaluateOperator<Eigen::SparseMatrix<double>>(ham_par.core(), ham_par.twoElectron(), container, diagonal_values);
    container.addToMatrix();
    return container.get_matrix();
}


/**
 *  Calculates sigma(pq) + sigma(qp)'s: all one-electron couplings for each annihilation-creation pair in the Fock space
 *  and stores them in sparse matrices for each pair combination
 *
 *  @return vector of sparse matrices containing the one-electron couplings for the (spin) Fock space
 *      Ordered as: sigma(00), sigma(01) + sigma(10), sigma(02)+ sigma(20), ...
 */
std::vector<Eigen::SparseMatrix<double>> FockSpace::calculateOneElectronCouplings() const {

    size_t K = this->get_K();
    size_t N = this->get_N();
    size_t dim = this->get_dimension();

    std::vector<std::vector<Eigen::Triplet<double>>> sparse_entries (K*(K+1)/2);
    std::vector<Eigen::SparseMatrix<double>> sparse_matrices (K*(K+1)/2, Eigen::SparseMatrix<double>(dim, dim));

    if (N == 0) {
        return sparse_matrices;
    }

    // Reserve appropriate amount of entries
    size_t reservation_size = FockSpace::calculateDimension(K-1, N-1);
    for (size_t p = 0; p < K; p++) {
        sparse_entries[p*(K+K+1-p)/2].reserve(reservation_size);
        for (size_t q = p + 1; q < K; q++) {
            sparse_entries[p*(K+K+1-p)/2 + q - p].reserve(2*reservation_size);
        }
    }

    ONV onv = this->makeONV(0);  // onv with address 0
    for (size_t I = 0; I < dim; I++) {  // I loops over all the addresses of the onv
        for (size_t e1 = 0; e1 < N; e1++) {  // e1 (electron 1) loops over the (number of) electrons
            size_t p = onv.get_occupation_index(e1);  // retrieve the index of a given electron
            // remove the weight from the initial address I, because we annihilate
            size_t address = I - this->get_vertex_weights(p, e1 + 1);
            // The e2 iteration counts the amount of encountered electrons for the creation operator
            // We only consider greater addresses than the initial one (because of symmetry)
            // Hence we only count electron after the annihilated electron (e1)
            sparse_entries[p*(K+K+1-p)/2].emplace_back(I, I, 1);
            size_t e2 = e1 + 1;
            size_t q = p + 1;
            int sign_e2 = 1;
            // perform a shift
            this->shiftUntilNextUnoccupiedOrbital<1>(onv, address, q, e2, sign_e2);
            size_t dex = 0;
            while (q < K) {
                size_t J = address + this->get_vertex_weights(q, e2);

                sparse_entries[p*(K+K+1-p)/2 + q - p].emplace_back(I, J, sign_e2);
                sparse_entries[p*(K+K+1-p)/2 + q - p].emplace_back(J, I, sign_e2);

                q++; // go to the next orbital
                dex++;
                // perform a shift
                this->shiftUntilNextUnoccupiedOrbital<1>(onv, address, q, e2, sign_e2);
            }  //  (creation)
        } // e1 loop (annihilation)

        // Prevent last permutation
        if (I < dim - 1) {
            this->setNextONV(onv);
        }
    }

    for (size_t k = 0; k < K*(K+1)/2 ; k++){
        sparse_matrices[k].setFromTriplets(sparse_entries[k].begin(), sparse_entries[k].end());
    }

    return sparse_matrices;
}



/**
 *  Evaluate the diagonal of the operator
 *
 *  @param one_op               the one-electron operator in an orthonormal orbital basis to be evaluated in the Fock space
 *
 *  @return the operator's diagonal evaluation in a vector with the dimension of the Fock space
 */
VectorX<double> FockSpace::evaluateOperatorDiagonal(const ScalarSQOneElectronOperator<double>& one_op) const {

    const auto K = one_op.dimension();

    if (K != this->K) {
        throw std::invalid_argument("FockSpace::evaluateOperatorDiagonal(ScalarSQOneElectronOperator<double>): Basis functions of the Fock space and the operator are incompatible.");
    }

    const auto& one_op_par = one_op.parameters();

    size_t N = this->get_N();
    size_t dim = this->get_dimension();

    VectorX<double> diagonal = VectorX<double>::Zero(dim);

    ONV onv = this->makeONV(0);  // onv with address 0
    for (size_t I = 0; I < dim; I++) {  // I loops over all addresses in the Fock space

        if (I > 0) {
            this->setNextONV(onv);
        }

        for (size_t e1 = 0; e1 < N; e1++) {  // A1 (annihilation 1)
            size_t p = onv.get_occupation_index(e1);
            diagonal(I) += one_op_par(p, p);
        }
    }

    return diagonal;
};


/**
 *  Evaluate the diagonal of the operator
 *
 *  @param two_op               the two-electron operator in an orthonormal orbital basis to be evaluated in the Fock space
 *
 *  @return the operator's diagonal evaluation in a vector with the dimension of the Fock space
 */
VectorX<double> FockSpace::evaluateOperatorDiagonal(const ScalarSQTwoElectronOperator<double>& two_op) const {

    auto K = two_op.dimension();
    if (K != this->K) {
        throw std::invalid_argument("FockSpace::evaluateOperatorDiagonal(ScalarSQTwoElectronOperator<double>): Basis functions of the Fock space and the operator are incompatible.");
    }

    size_t N = this->get_N();
    size_t dim = this->get_dimension();

    VectorX<double> diagonal = VectorX<double>::Zero(dim);

    const auto k_par = two_op.effectiveOneElectronPartition().parameters();
    const auto& two_op_par = two_op.parameters();

    ONV onv = this->makeONV(0);  // onv with address 0
    for (size_t I = 0; I < dim; I++) {  // I loops over all addresses in the Fock space

        if (I > 0) {
            this->setNextONV(onv);
        }

        for (size_t e1 = 0; e1 < N; e1++) {  // A1 (annihilation 1)
            size_t p = onv.get_occupation_index(e1);
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
 *  @param ham_par              the Hamiltonian parameters in an orthonormal orbital basis to be evaluated in the Fock space
 *
 *  @return the Hamiltonian's diagonal evaluation in a vector with the dimension of the Fock space
 */
VectorX<double> FockSpace::evaluateOperatorDiagonal(const SQHamiltonian<double>& ham_par) const {

    auto K = ham_par.dimension();
    if (K != this->K) {
        throw std::invalid_argument("FockSpace::evaluateOperatorDiagonal(SQHamiltonian<double>): Basis functions of the Fock space and the operator are incompatible.");
    }

    return this->evaluateOperatorDiagonal(ham_par.core()) + this->evaluateOperatorDiagonal(ham_par.twoElectron());
};



}  // namespace GQCP
