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
#include "QCModel/Geminals/APIGGeminalCoefficients.hpp"

#include "Mathematical/Representation/SquareMatrix.hpp"
#include "ONVBasis/SpinUnresolvedONVBasis.hpp"
#include "Utilities/miscellaneous.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param G                the APIG geminal coefficients
 */
APIGGeminalCoefficients::APIGGeminalCoefficients(const MatrixX<double>& G) :
    K {static_cast<size_t>(G.cols())},
    N_P {static_cast<size_t>(G.rows())},
    G {G} {}


/**
 *  Constructor that sets the geminal coefficients to zero
 *
 *  @param N_P      the number of electron pairs (= the number of geminals)
 *  @param K        the number of spatial orbitals
 */
APIGGeminalCoefficients::APIGGeminalCoefficients(const size_t N_P, const size_t K) :
    APIGGeminalCoefficients(MatrixX<double>::Zero(N_P, K)) {}


/*
 *  DESTRUCTOR
 */

APIGGeminalCoefficients::~APIGGeminalCoefficients() {}


/*
 *  OPERATORS
 */

/**
 *  @param i            the zero-based index of the geminal, i.e. subscript of the geminal coefficient: i is in [0, N_P[ with N_P the number of electron pairs
 *  @param p            the zero-based index of the orbital, i.e. superscript of the geminal coefficient
 * 
 *  @return an element of the AP1roG geminal coefficient matrix G_i^a
 */
double APIGGeminalCoefficients::operator()(const size_t i, const size_t p) const {

    return this->G(i, p);
}


/*
 *  NAMED CONSTRUCTORS
 */

/**
 *  @param g        the geminal coefficients in a vector representation that is in column-major storage
 *  @param N_P      the number of electron pairs (= the number of geminals)
 *  @param K        the number of spatial orbitals
 */
APIGGeminalCoefficients APIGGeminalCoefficients::FromColumnMajor(const VectorX<double>& g, const size_t N_P, const size_t K) {

    return APIGGeminalCoefficients(MatrixX<double>::FromColumnMajorVector(g, N_P, K));
}


/**
 *  @param g        the geminal coefficients in a vector representation that is in row-major storage
 *  @param N_P      the number of electron pairs (= the number of geminals)
 *  @param K        the number of spatial orbitals
 */
APIGGeminalCoefficients APIGGeminalCoefficients::FromRowMajor(const VectorX<double>& g, const size_t N_P, const size_t K) {

    return APIGGeminalCoefficients(MatrixX<double>::FromRowMajorVector(g, N_P, K));
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
size_t APIGGeminalCoefficients::numberOfGeminalCoefficients(size_t N_P, size_t K) {

    // Check if we can have N_P geminals in K orbitals
    if (N_P >= K) {
        throw std::invalid_argument("APIGGeminalCoefficients::numberOfGeminalCoefficients(size_t, size_t): Can't have that many geminals in this few number of orbitals.");
    }

    return N_P * K;
}


/*
 *  PUBLIC METHODS
 */

/**
 *  @param onv      the doubly-occupied (spin-resolved) ONV that is being projected on
 *
 *  @return the overlap of the APIG wave function with the given ONV, i.e. the projection of the APIG wave function onto that ONV
 */
double APIGGeminalCoefficients::overlap(const SpinUnresolvedONV& onv) const {

    // For a pure APIG, the formula has to be the most general one


    // Construct the matrix G(m) which only has the occupied columns of G in the given doubly-occupied (spin-resolved) ONV 'm'
    MatrixX<double> G = this->asMatrix();  // geminal coefficients as a matrix
    SquareMatrix<double> Gm = SquareMatrix<double>::Zero(this->N_P, this->N_P);

    // TODO: wait until the syntax G(Eigen::placeholders::all, occupation_indices) is released in a stable Eigen release
    for (size_t e = 0; e < this->N_P; e++) {  // loop over all electrons
        size_t occupation_index = onv.get_occupation_index(e);

        Gm.col(e) = G.col(occupation_index);
    }


    // Calculate the permanent of Gm to obtain the coefficient
    return Gm.permanent_ryser();
}


}  // namespace GQCP
