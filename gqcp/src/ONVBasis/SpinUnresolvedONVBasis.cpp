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
 *  MARK: Constructors
 */

/**
 *  @param M            The number of spinors/spin-orbitals.
 *  @param N            The number of electrons, i.e. the number of occupied spinors/spin-orbitals.
 */
SpinUnresolvedONVBasis::SpinUnresolvedONVBasis(const size_t M, const size_t N) :
    M {M},
    N {N} {

    // Set up the vertex weights for the addressing scheme for a full spin-unresolved ONV basis. This addressing scheme is taken from Helgaker, Jørgensen, Olsen (2000).

    // Create a zero matrix of dimensions (M+1)x(N+1)
    this->vertex_weights = std::vector<std::vector<size_t>>(M + 1, std::vector<size_t>(N + 1, 0));

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
    for (size_t p = 0; p < M - N + 1; p++) {
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

    for (size_t m = 1; m < N + 1; m++) {
        for (size_t p = m; p < (M - N + m) + 1; p++) {
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
size_t SpinUnresolvedONVBasis::calculateDimension(const size_t M, const size_t N) {
    const auto dim_double = boost::math::binomial_coefficient<double>(static_cast<unsigned>(M), static_cast<unsigned>(N));
    try {
        return boost::numeric::converter<size_t, double>::convert(dim_double);
    } catch (boost::numeric::bad_numeric_cast& e) {
        throw std::overflow_error("SpinUnresolvedONVBasis::calculateDimension(const size_t, const size_t): " + std::string(e.what()));
    }
}


/*
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
size_t SpinUnresolvedONVBasis::arcWeight(const size_t p, const size_t n) const {

    // Arc weights and vertex weights are related. This relation is found in Helgaker, Jørgensen, Olsen (2000) in chapter 11.3.6.
    return this->vertexWeight(p, n + 1);
}


/**
 *  Calculate the address (i.e. the ordering number) of an unsigned representation of a spin-unresolved ONV.
 * 
 *  @param representation      The unsigned representation of a spin-unresolved ONV.
 *
 *  @return The address corresponding to the unsigned representation of a spin-unresolved ONV.
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
size_t SpinUnresolvedONVBasis::nextPermutationOf(const size_t representation) const {

    // t gets this->representation's least significant 0 bits set to 1
    unsigned long t = representation | (representation - 1UL);

    // Next set to 1 the most significant bit to change,
    // set to 0 the least significant ones, and add the necessary 1 bits.
    return (t + 1UL) | (((~t & (t + 1UL)) - 1UL) >> (__builtin_ctzl(representation) + 1UL));
}


/**
 *  Calculate the unsigned representation of a spin-unresolved ONV that corresponds to the address/ordering number in this ONV basis.
 *
 *  @param address                 The address/ordering number of a spin-unresolved ONV in this ONV basis.
 *
 *  @return The unsigned representation of a spin-unresolved ONV that corresponds to the address/ordering number in this ONV basis.
 */
size_t SpinUnresolvedONVBasis::representationOf(size_t address) const {
    size_t representation = 0;
    if (this->numberOfElectrons() != 0) {
        representation = 0;
        size_t m = this->numberOfElectrons();  // counts the number of electrons in the spin string up to orbital p

        for (size_t p = this->numberOfOrbitals(); p > 0; p--) {  // p is an orbital index
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
 *  Create the ONV that corresponds to the given address in this ONV basis.
 * 
 *  @param address                 The address/ordering number of a spin-unresolved ONV in this ONV basis.
 *
 *  @return The ONV that corresponds to the given address in this ONV basis.
 */
SpinUnresolvedONV SpinUnresolvedONVBasis::constructONVFromAddress(const size_t address) const {

    SpinUnresolvedONV onv {this->numberOfOrbitals(), this->numberOfElectrons()};
    this->transformONVCorrespondingToAddress(onv, address);
    return onv;
}


/**
 *  Modify a `SpinResolvedONV` to the next allowed ONV in this ONV basis.
 * 
 *  @param onv      A spin-resolved ONV.
 */
void SpinUnresolvedONVBasis::transformONVToNextPermutation(SpinUnresolvedONV& onv) const {

    onv.replaceRepresentationWith(this->nextPermutationOf(onv.unsignedRepresentation()));
}


/**
 *  Modify a `SpinResolvedONV` to the one with the given address in this ONV basis.
 *
 *  @param onv          A spin-resolved ONV.
 *  @param address      The target address in this ONV basis.
 */
void SpinUnresolvedONVBasis::transformONVCorrespondingToAddress(SpinUnresolvedONV& onv, const size_t address) const {

    onv.replaceRepresentationWith(this->representationOf(address));
}


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
size_t SpinUnresolvedONVBasis::countOneElectronCouplings(const SpinUnresolvedONV& onv) const {

    size_t V = this->numberOfOrbitals() - this->numberOfElectrons();  // number of virtual orbitals
    size_t coupling_count = 0;

    for (size_t e1 = 0; e1 < this->numberOfElectrons(); e1++) {
        size_t p = onv.occupationIndexOf(e1);
        coupling_count += (V + e1 - p);  // number of virtuals with an index larger than p
    }

    return coupling_count;
}


/**
 *  @return The total number of non-zero and non-diagonal couplings of a one-electron operator in this ONV basis.
 */
size_t SpinUnresolvedONVBasis::countTotalOneElectronCouplings() const {

    const auto M = this->numberOfOrbitals();
    const auto N = this->numberOfElectrons();

    return (M - N) * N * (this->dimension());
}


/**
 *  @return The total number of non-zero and non-diagonal couplings of a two-electron operator in this ONV basis.
 */
size_t SpinUnresolvedONVBasis::countTotalTwoElectronCouplings() const {

    const auto M = this->numberOfOrbitals();
    const auto N = this->numberOfElectrons();

    size_t two_electron_permutation = 0;  // all distributions for two electrons over the virtual orbitals
    if (M - N >= 2) {
        two_electron_permutation = SpinUnresolvedONVBasis::calculateDimension(M - N, 2) * N * (N - 1) * (this->dimension()) / 2;
    }

    return two_electron_permutation + countTotalOneElectronCouplings();
}


/**
 *  Calculate the number of ONVs (with a larger address) that a given spin-unresolved ONV would couple with in this ONV basis, through a two-electron operator.
 * 
 *  @param onv          The spin-unresolved ONV.
 *
 *  @return The number of ONVs (with a larger address) that a given spin-unresolved ONV would couple with in this ONV basis, through a two-electron operator.
 */
size_t SpinUnresolvedONVBasis::countTwoElectronCouplings(const SpinUnresolvedONV& onv) const {

    const auto M = this->numberOfOrbitals();
    const auto N = this->numberOfElectrons();

    size_t V = M - N;  // number of virtual orbitals
    size_t coupling_count = 0;

    for (size_t e1 = 0; e1 < N; e1++) {

        size_t p = onv.occupationIndexOf(e1);
        coupling_count += (V + e1 - p);  //  one electron part

        for (size_t e2 = e1 + 1; e2 < N; e2++) {

            size_t q = onv.occupationIndexOf(e2);
            size_t coupling_count2 = (V + e2 - q);
            coupling_count += (V - coupling_count2) * coupling_count2;

            if (coupling_count2 > 1) {
                coupling_count += SpinUnresolvedONVBasis::calculateDimension(coupling_count2, 2);
            }
        }
    }

    return coupling_count;
}


/**
 *  Calculate all one-electron coupling elements for this spin-unresolved ONV basis. These are all the intermediate matrices sigma(pq)_{IJ}, as defined by Helgaker, Jørgensen, Olsen (2000).
 *
 *  @return A vector of sparse matrices containing the one-electron coupling elements for this spin-unresolved ONV basis. The elements of this vector are ordered through the one-electron excitation (pq) in ascending order: sigma(00), sigma(01) + sigma(10), sigma(02)+ sigma(20), ...
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
            while (q < K) {
                size_t J = address + this->vertexWeight(q, e2);

                sparse_entries[p * (K + K + 1 - p) / 2 + q - p].emplace_back(I, J, sign_e2);
                sparse_entries[p * (K + K + 1 - p) / 2 + q - p].emplace_back(J, I, sign_e2);

                q++;  // go to the next orbital
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


/*
 *  MARK: Iterating
 */

/**
 *  Iterate over all ONVs in this ONV basis and apply the given callback function.
 * 
 *  @param callback            The function to be applied in every iteration. Its supplied arguments are a spin-unresolved ONV and its corresponding addresses.
 */
void SpinUnresolvedONVBasis::forEach(const std::function<void(const SpinUnresolvedONV&, const size_t)>& callback) const {

    const auto dim = this->dimension();
    SpinUnresolvedONV onv = this->constructONVFromAddress(0);

    // Loop over all addresses of the ONVs.
    for (size_t I = 0; I < dim; I++) {

        callback(onv, I);

        // Prevent the last permutation from occurring, as this would cause errors.
        if (I < dim - 1) {
            this->transformONVToNextPermutation(onv);
        }
    }
}


/*
 *  MARK: Dense generalized operator evaluations
 */

/**
 *  Calculate the dense matrix representation of a generalized two-electron operator in this ONV basis.
 *
 *  @param g                A generalized two-electron operator expressed in an orthonormal orbital basis.
 *
 *  @return A dense matrix represention of the two-electron operator.
 */
SquareMatrix<double> SpinUnresolvedONVBasis::evaluateOperatorDense(const ScalarGSQTwoElectronOperator<double>& g) const {

    // In order to avoid duplicate code, we choose to delegate this method to the evaluation of a `GSQHamiltonian` that contains no core contributions. This does not affect performance significantly, because the bottleneck will always be the iteration over the whole ONV basis.
    const auto zero = ScalarGSQOneElectronOperator<double>::Zero(g.numberOfOrbitals());
    const GSQHamiltonian<double> hamiltonian {zero, g};

    return this->evaluateOperatorDense(hamiltonian);
}


/**
 *  Calculate the dense matrix representation of a generalized Hamiltonian in this ONV basis.
 *
 *  @param hamiltonian      A generalized Hamiltonian expressed in an orthonormal orbital basis.
 *
 *  @return A dense matrix represention of the Hamiltonian.
 */
SquareMatrix<double> SpinUnresolvedONVBasis::evaluateOperatorDense(const GSQHamiltonian<double>& hamiltonian) const {

    if (hamiltonian.numberOfOrbitals() != this->numberOfOrbitals()) {
        throw std::invalid_argument("SpinUnresolvedONVBasis::evaluateOperatorDense(const USQHamiltonian<double>&): The number of orbitals of this ONV basis and the given Hamiltonian are incompatible.");
    }

    // Initialize a container for the dense matrix representation, and fill it with the general evaluation function.
    MatrixRepresentationEvaluationContainer<SquareMatrix<double>> container {this->dimension()};
    this->evaluate<SquareMatrix<double>>(hamiltonian, container);

    return container.evaluation();
}


/*
 *  MARK: Dense unrestricted operator evaluations
 */

/**
 *  Calculate the dense matrix representation of a component of an unrestricted one-electron operator in this ONV basis.
 *
 *  @param f                A component of an unrestricted one-electron operator expressed in an orthonormal orbital basis.
 *
 *  @return A dense matrix represention of the one-electron operator.
 */
SquareMatrix<double> SpinUnresolvedONVBasis::evaluateOperatorDense(const ScalarUSQOneElectronOperatorComponent<double>& f) const {

    // We may convert an unrestricted component into the generalized representation.
    const auto f_generalized = ScalarGSQOneElectronOperator<double>::FromUnrestrictedComponent(f);
    return this->evaluateOperatorDense(f_generalized);
}


/**
 *  Calculate the dense matrix representation of a component of an unrestricted two-electron operator in this ONV basis.
 *
 *  @param g                A component of an unrestricted two-electron operator expressed in an orthonormal orbital basis.
 *
 *  @return A dense matrix represention of the one-electron operator.
 */
SquareMatrix<double> SpinUnresolvedONVBasis::evaluateOperatorDense(const ScalarPureUSQTwoElectronOperatorComponent<double>& g) const {

    // We may convert an unrestricted component into the generalized representation.
    const auto g_generalized = ScalarGSQTwoElectronOperator<double>::FromUnrestrictedComponent(g);
    return this->evaluateOperatorDense(g_generalized);
}


/*
 *  MARK: Diagonal generalized operator evaluations
 */

/**
 *  Calculate the diagonal of the matrix representation of a generalized one-electron operator in this ONV basis.
 *
 *  @param f_op             A generalized one-electron operator expressed in an orthonormal orbital basis.
 *
 *  @return The diagonal of the dense matrix represention of the one-electron operator.
 */
VectorX<double> SpinUnresolvedONVBasis::evaluateOperatorDiagonal(const ScalarGSQOneElectronOperator<double>& f_op) const {

    const auto K = f_op.numberOfOrbitals();

    if (K != this->numberOfOrbitals()) {
        throw std::invalid_argument("SpinUnresolvedONVBasis::evaluateOperatorDiagonal(const ScalarGSQOneElectronOperator<double>&): The number of orbitals of this ONV basis and the given operator are incompatible.");
    }


    // Prepare some variables.
    const auto& f = f_op.parameters();

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
            diagonal(I) += f(p, p);
        }
    }

    return diagonal;
};


/**
 *  Calculate the diagonal of the matrix representation of a generalized two-electron operator in this ONV basis.
 *
 *  @param g_op             A generalized two-electron operator expressed in an orthonormal orbital basis.
 *
 *  @return The diagonal of the dense matrix represention of the two-electron operator.
 */
VectorX<double> SpinUnresolvedONVBasis::evaluateOperatorDiagonal(const ScalarGSQTwoElectronOperator<double>& g_op) const {

    auto K = g_op.numberOfOrbitals();
    if (K != this->numberOfOrbitals()) {
        throw std::invalid_argument("SpinUnresolvedONVBasis::evaluateOperatorDiagonal(const ScalarGSQTwoElectronOperator<double>&): The number of orbitals of this ONV basis and the given operator are incompatible.");
    }

    size_t N = this->numberOfElectrons();
    size_t dim = this->dimension();

    VectorX<double> diagonal = VectorX<double>::Zero(dim);

    const auto k = g_op.effectiveOneElectronPartition().parameters();
    const auto& g = g_op.parameters();

    SpinUnresolvedONV onv = this->constructONVFromAddress(0);  // onv with address 0
    for (size_t I = 0; I < dim; I++) {                         // I loops over all addresses in this ONV basis

        if (I > 0) {
            this->transformONVToNextPermutation(onv);
        }

        for (size_t e1 = 0; e1 < N; e1++) {  // A1 (annihilation 1)
            size_t p = onv.occupationIndexOf(e1);
            diagonal(I) += k(p, p);

            for (size_t q = 0; q < K; q++) {  // q loops over SOs
                if (onv.isOccupied(q)) {
                    diagonal(I) += 0.5 * g(p, p, q, q);
                } else {
                    diagonal(I) += 0.5 * g(p, q, q, p);
                }
            }
        }
    }

    return diagonal;
};


/**
 *  Calculate the diagonal of the dense matrix representation of a generalized Hamiltonian in this ONV basis.
 *
 *  @param hamiltonian      A generalized Hamiltonian expressed in an orthonormal orbital basis.
 *
 *  @return The diagonal of the dense matrix represention of the Hamiltonian.
 */
VectorX<double> SpinUnresolvedONVBasis::evaluateOperatorDiagonal(const GSQHamiltonian<double>& sq_hamiltonian) const {

    return this->evaluateOperatorDiagonal(sq_hamiltonian.core()) + this->evaluateOperatorDiagonal(sq_hamiltonian.twoElectron());
};


/*
 *  MARK: Sparse generalized operator evaluations
 */

/**
 *  Calculate the sparse matrix representation of a generalized one-electron operator in this ONV basis.
 *
 *  @param f                A generalized one-electron operator expressed in an orthonormal orbital basis.
 *
 *  @return A sparse matrix represention of the one-electron operator.
 */
Eigen::SparseMatrix<double> SpinUnresolvedONVBasis::evaluateOperatorSparse(const ScalarGSQOneElectronOperator<double>& f) const {

    if (f.numberOfOrbitals() != this->numberOfOrbitals()) {
        throw std::invalid_argument("SpinUnresolvedONVBasis::evaluateOperatorSparse(const ScalarGSQOneElectronOperator<double>&): The number of orbitals of the ONV basis and the operator are incompatible.");
    }

    // Initialize a container for the sparse matrix representation, and reserve an appropriate amount of memory for it.
    MatrixRepresentationEvaluationContainer<Eigen::SparseMatrix<double>> container {this->dimension()};
    size_t memory = this->dimension() + this->countTotalOneElectronCouplings();
    container.reserve(memory);

    // Evaluate the one-electron operator and add the evaluations to the sparse matrix representation.
    this->evaluate<Eigen::SparseMatrix<double>>(f, container);

    // Finalize the creation of the sparse matrix and return the result.
    container.addToMatrix();
    return container.evaluation();
}


/**
 *  Calculate the sparse matrix representation of a generalized two-electron operator in this ONV basis.
 *
 *  @param g                A generalized two-electron operator expressed in an orthonormal orbital basis.
 *
 *  @return A sparse matrix represention of the two-electron operator.
 */
Eigen::SparseMatrix<double> SpinUnresolvedONVBasis::evaluateOperatorSparse(const ScalarGSQTwoElectronOperator<double>& g) const {

    // In order to avoid duplicate code, we choose to delegate this method to the evaluation of a `GSQHamiltonian` that contains no core contributions. This does not affect performance significantly, because the bottleneck will always be the iteration over the whole ONV basis.
    const auto zero = ScalarGSQOneElectronOperator<double>::Zero(g.numberOfOrbitals());
    const GSQHamiltonian<double> hamiltonian {zero, g};

    return this->evaluateOperatorSparse(hamiltonian);
}


/**
 *  Calculate the sparse matrix representation of a generalized Hamiltonian in this ONV basis.
 *
 *  @param hamiltonian      A generalized Hamiltonian expressed in an orthonormal orbital basis.
 *
 *  @return A sparse matrix represention of the Hamiltonian.
 */
Eigen::SparseMatrix<double> SpinUnresolvedONVBasis::evaluateOperatorSparse(const GSQHamiltonian<double>& hamiltonian) const {

    if (hamiltonian.numberOfOrbitals() != this->numberOfOrbitals()) {
        throw std::invalid_argument("SpinUnresolvedONVBasis::evaluateOperatorSparse(const GSQHamiltonian<double>& hamiltonian): The number of orbitals of this ONV basis and the given Hamiltonian are incompatible.");
    }

    // Initialize a container for the sparse matrix representation, and reserve an appropriate amount of memory for it.
    MatrixRepresentationEvaluationContainer<Eigen::SparseMatrix<double>> container {this->dimension()};
    size_t memory = this->dimension() + this->countTotalTwoElectronCouplings();
    container.reserve(memory);

    // Evaluate the Hamiltonian and add the evaluations to the sparse matrix representation.
    this->evaluate<Eigen::SparseMatrix<double>>(hamiltonian, container);

    // Finalize the creation of the sparse matrix and return the result.
    container.addToMatrix();
    return container.evaluation();
}


/*
 *  MARK: Sparse unrestricted operator evaluations
 */

/**
 *  Calculate the sparse matrix representation of a component of an unrestricted one-electron operator in this ONV basis.
 *
 *  @param f                A component of an unrestricted one-electron operator expressed in an orthonormal orbital basis.
 *
 *  @return A sparse matrix represention of the one-electron operator.
 */
Eigen::SparseMatrix<double> SpinUnresolvedONVBasis::evaluateOperatorSparse(const ScalarUSQOneElectronOperatorComponent<double>& f) const {

    // We may convert an unrestricted component into the generalized representation.
    const auto f_generalized = ScalarGSQOneElectronOperator<double>::FromUnrestrictedComponent(f);
    return this->evaluateOperatorSparse(f_generalized);
}


/**
 *  Calculate the sparse matrix representation of a component of an unrestricted two-electron operator in this ONV basis.
 *
 *  @param g                A component of an unrestricted two-electron operator expressed in an orthonormal orbital basis.
 *
 *  @return A sparse matrix represention of the one-electron operator.
 */
Eigen::SparseMatrix<double> SpinUnresolvedONVBasis::evaluateOperatorSparse(const ScalarPureUSQTwoElectronOperatorComponent<double>& g) const {

    // We may convert an unrestricted component into the generalized representation.
    const auto g_generalized = ScalarGSQTwoElectronOperator<double>::FromUnrestrictedComponent(g);
    return this->evaluateOperatorSparse(g_generalized);
}


/*
 *  MARK: Generalized matrix-vector product evaluations
 */

/**
 *  Calculate the matrix-vector product of (the matrix representation of) a generalized one-electron operator with the given coefficient vector.
 *
 *  @param f                A generalized one-electron operator expressed in an orthonormal orbital basis.
 *  @param x                The coefficient vector of a linear expansion.
 *
 *  @return The coefficient vector of the linear expansion after being acted on with the given (matrix representation of) the one-electron operator.
 */
VectorX<double> SpinUnresolvedONVBasis::evaluateOperatorMatrixVectorProduct(const ScalarGSQOneElectronOperator<double>& f, const VectorX<double>& x) const {

    if (f.numberOfOrbitals() != this->numberOfOrbitals()) {
        throw std::invalid_argument("SpinUnresolvedONVBasis::evaluateOperatorMatrixVectorProduct(const ScalarGSQOneElectronOperator<double>&, const VectorX<double>&): The number of orbitals of this ONV basis and the operator are incompatible.");
    }

    // Initialize a container for the matrix-vector product, and fill it with the general evaluation function.
    MatrixRepresentationEvaluationContainer<VectorX<double>> container {x};
    this->evaluate<VectorX<double>>(f, container);

    return container.evaluation();
}


/**
 *  Calculate the matrix-vector product of (the matrix representation of) a generalized two-electron operator with the given coefficient vector.
 *
 *  @param g                A generalized two-electron operator expressed in an orthonormal orbital basis.
 *  @param x                The coefficient vector of a linear expansion.
 *
 *  @return The coefficient vector of the linear expansion after being acted on with the given (matrix representation of) the two-electron operator.
 */
VectorX<double> SpinUnresolvedONVBasis::evaluateOperatorMatrixVectorProduct(const ScalarGSQTwoElectronOperator<double>& g, const VectorX<double>& x) const {

    // In order to avoid duplicate code, we choose to delegate this method to the evaluation of a `GSQHamiltonian` that contains no core contributions. This does not affect performance significantly, because the bottleneck will always be the iteration over the whole ONV basis.
    const auto zero = ScalarGSQOneElectronOperator<double>::Zero(g.numberOfOrbitals());
    const GSQHamiltonian<double> hamiltonian {zero, g};

    return this->evaluateOperatorMatrixVectorProduct(hamiltonian, x);
}


/**
 *  Calculate the matrix-vector product of (the matrix representation of) a generalized Hamiltonian with the given coefficient vector.
 *
 *  @param hamiltonian      A generalized Hamiltonian expressed in an orthonormal orbital basis.
 *  @param x                The coefficient vector of a linear expansion.
 *
 *  @return The coefficient vector of the linear expansion after being acted on with the given (matrix representation of) the Hamiltonian.
 */
VectorX<double> SpinUnresolvedONVBasis::evaluateOperatorMatrixVectorProduct(const GSQHamiltonian<double>& hamiltonian, const VectorX<double>& x) const {

    if (hamiltonian.numberOfOrbitals() != this->numberOfOrbitals()) {
        throw std::invalid_argument("SpinUnresolvedONVBasis::evaluateOperatorMatrixVectorProduct(const USQHamiltonian<double>&, const VectorX<double>& x): The number of orbitals of this ONV basis and the given Hamiltonian are incompatible.");
    }

    // Initialize a container for the matrix-vector product, and fill it with the general evaluation function.
    MatrixRepresentationEvaluationContainer<VectorX<double>> container {x};
    this->evaluate<VectorX<double>>(hamiltonian, container);

    return container.evaluation();
}


}  // namespace GQCP
