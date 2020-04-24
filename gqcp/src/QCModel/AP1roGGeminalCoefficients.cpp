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

#include "QCModel/Geminals/AP1roGGeminalCoefficients.hpp"

#include "ONVBasis/SpinUnresolvedONVBasis.hpp"
#include "QCModel/Geminals/APIGGeminalCoefficients.hpp"
#include "Utilities/miscellaneous.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */


/**
 *  @param G            the AP1roG geminal coefficients (not including the identity matrix on the left), as a block matrix
 */
AP1roGGeminalCoefficients::AP1roGGeminalCoefficients(const BlockMatrix<double>& G, const size_t N_P, const size_t K) :
    N_P {N_P},
    K {K},
    G {G} {}


/**
 *  @param G            the AP1roG geminal coefficients (not including the identity matrix on the left)
 */
AP1roGGeminalCoefficients::AP1roGGeminalCoefficients(const MatrixX<double>& G) :
    // Don't use a delegated constructor for code readability
    N_P {static_cast<size_t>(G.rows())},
    K {static_cast<size_t>(G.rows() + G.cols())},
    G {BlockMatrix<double>(0, this->N_P, this->N_P, this->K, G)} {}


/**
 *  Constructor that sets the geminal coefficients to zero
 *
 *  @param N_P      the number of electron pairs (= the number of geminals)
 *  @param K        the number of spatial orbitals
 */
AP1roGGeminalCoefficients::AP1roGGeminalCoefficients(const size_t N_P, const size_t K) :
    AP1roGGeminalCoefficients(MatrixX<double>::Zero(N_P, K - N_P)) {}


/*
 *  DESTRUCTOR
 */

AP1roGGeminalCoefficients::~AP1roGGeminalCoefficients() {}


/*
 *  OPERATORS
 */

/**
 *  @param i            the zero-based index of the geminal, i.e. subscript of the geminal coefficient: i is in [0, N_P[ with N_P the number of electron pairs
 *  @param a            the zero-based index of the occupied orbital, i.e. superscript of the geminal coefficient: a is in [N_P, K[ with K the number of spatial orbitals
 * 
 *  @return an element of the AP1roG geminal coefficient matrix G_i^a
 */
double AP1roGGeminalCoefficients::operator()(const size_t i, const size_t a) const {

    return this->G(i, a);  // BlockMatrix implements operator() as we would expect
}


/*
 *  NAMED CONSTRUCTORS
 */

/**
 *  @param sq_hamiltonian       the Hamiltonian expressed in an orthonormal basis
 *  @param N_P                  the number of orbitals
 *
 *  @return the AP1roG geminal coefficients in the weak interaction limit
 */
AP1roGGeminalCoefficients AP1roGGeminalCoefficients::WeakInteractionLimit(const SQHamiltonian<double>& sq_hamiltonian, const size_t N_P) {

    const auto K = sq_hamiltonian.dimension();
    const auto number_of_geminal_coefficients = AP1roGGeminalCoefficients::numberOfGeminalCoefficients(N_P, K);
    const auto& h = sq_hamiltonian.core().parameters();         // core Hamiltonian integrals
    const auto& g = sq_hamiltonian.twoElectron().parameters();  // two-electron integrals

    // Provide the weak interaction limit values for the geminal coefficients
    BlockMatrix<double> G {0, N_P, N_P, K};
    for (size_t i = 0; i < N_P; i++) {
        for (size_t a = N_P; a < K; a++) {
            G(i, a) = -g(a, i, a, i) / (2 * (h(a, a) - h(i, i)));
        }
    }

    return AP1roGGeminalCoefficients(G, N_P, K);
}


/**
 *  @param g        the geminal coefficients in a vector representation that is in column-major storage
 *
 *  @param N_P      the number of electron pairs (= the number of geminals)
 *  @param K        the number of spatial orbitals
 */
AP1roGGeminalCoefficients AP1roGGeminalCoefficients::FromColumnMajor(const VectorX<double>& g, size_t N_P, size_t K) {

    const size_t rows = N_P;
    const size_t cols = K - N_P;  // the total number of columns in the (reduced) AP1roG geminal coefficient matrix

    const MatrixX<double> M = MatrixX<double>::FromColumnMajorVector(g, rows, cols);  // the block of the actual entries of the geminal coefficient matrix

    return AP1roGGeminalCoefficients(BlockMatrix<double>(0, N_P, N_P, K, M),
                                     N_P, K);  // an encapsulating object that implements operator() in an intuitive way
}


/**
 *  @param g        the geminal coefficients in a vector representation that is in row-major storage
 *
 *  @param N_P      the number of electron pairs (= the number of geminals)
 *  @param K        the number of spatial orbitals
 */
AP1roGGeminalCoefficients AP1roGGeminalCoefficients::FromRowMajor(const VectorX<double>& g, const size_t N_P, const size_t K) {

    const size_t rows = N_P;
    const size_t cols = K - N_P;  // the total number of columns in the (reduced) AP1roG geminal coefficient matrix

    const MatrixX<double> M = MatrixX<double>::FromRowMajorVector(g, rows, cols);  // the block of the actual entries of the geminal coefficient matrix

    return AP1roGGeminalCoefficients(BlockMatrix<double>(0, N_P, N_P, K, M), N_P, K);  // an encapsulating object that implements operator() in an intuitive way
}


/*
 *  STATIC PUBLIC METHODS
 */

/**
 *  @param N_P      the number of electron pairs (= the number of geminals)
 *  @param K        the number of spatial orbitals
 *
 *  @return the number of 'free' geminal coefficients
 */
size_t AP1roGGeminalCoefficients::numberOfGeminalCoefficients(const size_t N_P, const size_t K) {

    // Check if we can have N_P geminals in K orbitals
    if (N_P >= K) {
        throw std::invalid_argument("AP1roGGeminalCoefficients::numberOfVariables(size_t, size_t): Can't have that many geminals in this few number of orbitals.");
    }

    return N_P * (K - N_P);
}


/*
 *  PUBLIC METHODS
 */

/**
 *  @return the total geminal coefficient matrix, including the identity matrix block
 */
MatrixX<double> AP1roGGeminalCoefficients::asMatrix() const {

    // Initialize the total geminal coefficient matrix
    MatrixX<double> G_total = MatrixX<double>::Zero(this->N_P, this->K);

    // The AP1roG coefficients are the identity matrix in the leftmost (N_P x N_P)-block
    G_total.topLeftCorner(this->N_P, this->N_P) = MatrixX<double>::Identity(this->N_P, this->N_P);

    // Set the right AP1roG coefficient block
    G_total.topRightCorner(this->N_P, this->K - this->N_P) = this->G.asMatrix();

    return G_total;
}


/**
 *  @return the geminal coefficients as a row-major vector, excluding the identity block
 */
VectorX<double> AP1roGGeminalCoefficients::asVector() const {

    return this->G.asVector();
}


/**
 *  @param onv      the doubly-occupied (spin-resolved) ONV that is being projected on
 *
 *  @return the overlap of the AP1roG wave function with the given ONV, i.e. the projection of the APIG wave function onto that ONV
 */
double AP1roGGeminalCoefficients::overlap(const SpinUnresolvedONV& onv) const {

    // For an AP1roG wave function, we use a simplification for singly and doubly pair-excited ONVs

    SpinUnresolvedONVBasis onv_basis {this->K, this->N_P};  // the doubly-occupied spin-resolved ONV basis
    SpinUnresolvedONV reference = onv_basis.makeONV(0);

    if (onv.countNumberOfDifferences(reference) == 0) {  // no excitations
        return 1.0;
    }

    else if (onv.countNumberOfDifferences(reference) == 2) {  // one pair excitation

        size_t i = reference.findDifferentOccupations(onv)[0];
        size_t a = onv.findDifferentOccupations(reference)[0];

        return this->operator()(i, a);
    }

    else if (onv.countNumberOfDifferences(reference) == 4) {  // two pair excitations

        auto different_occupied = reference.findDifferentOccupations(onv);
        auto different_virtual = onv.findDifferentOccupations(reference);

        size_t i = different_occupied[0];
        size_t j = different_occupied[1];
        size_t a = different_virtual[0];
        size_t b = different_virtual[1];

        return this->operator()(i, a) * this->operator()(j, b) + this->operator()(j, a) * this->operator()(i, b);
    }

    else {  // use the general formula if the difference is more than two pair excitations

        APIGGeminalCoefficients APIG {this->asMatrix()};
        return APIG.overlap(onv);
    }
}


}  // namespace GQCP
