// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2018  the GQCG developers
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
#include "geminals/APIGGeminalCoefficients.hpp"

#include "FockSpace/FockSpace.hpp"
#include "utilities/miscellaneous.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  Default constructor setting everything to zero
 */
APIGGeminalCoefficients::APIGGeminalCoefficients() :
    BaseAPIGGeminalCoefficients()
{}


/**
 *  @param g        the geminal coefficients in a vector representation that is in row-major storage
 *
 *  @param N_P      the number of electron pairs (= the number of geminals)
 *  @param K        the number of spatial orbitals
 */
APIGGeminalCoefficients::APIGGeminalCoefficients(const Eigen::VectorXd& g, size_t N_P, size_t K) :
    BaseAPIGGeminalCoefficients(g, N_P, K)
{
    if (APIGGeminalCoefficients::numberOfGeminalCoefficients(N_P, K) != g.size()) {
        throw std::invalid_argument("The specified N_P and K are not compatible with the given vector of geminal coefficients.");
    }
}


/**
 *  Constructor that sets the geminal coefficients to zero
 *
 *  @param N_P      the number of electron pairs (= the number of geminals)
 *  @param K        the number of spatial orbitals
 */
APIGGeminalCoefficients::APIGGeminalCoefficients(size_t N_P, size_t K) :
    APIGGeminalCoefficients(Eigen::VectorXd::Zero(APIGGeminalCoefficients::numberOfGeminalCoefficients(N_P, K)), N_P, K)
{}


/**
 *  @param G        the geminal coefficients in a matrix representation
 */
APIGGeminalCoefficients::APIGGeminalCoefficients(const Eigen::MatrixXd& G) :
    BaseAPIGGeminalCoefficients()
{

    Eigen::MatrixXd G_transpose = G.transpose();

    this->g = Eigen::Map<const Eigen::VectorXd>(G_transpose.data(), G_transpose.cols()*G_transpose.rows());
    this->K = G.cols();
    this->N_P = G.rows();
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
        throw std::invalid_argument("Can't have that many geminals in this few number of orbitals.");
    }

    return N_P * K;
}



/*
 *  PUBLIC METHODS
 */

/**
 *  @return the geminal coefficients in matrix form
 */
Eigen::MatrixXd APIGGeminalCoefficients::asMatrix() const {

    using RowMajorMatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    Eigen::RowVectorXd g_row = this->g;

    return Eigen::Map<RowMajorMatrixXd, Eigen::RowMajor>(g_row.data(), this->N_P, this->K);
}


/**
 *  @param vector_index     the vector index of the geminal coefficient
 *
 *  @return the major (geminal, subscript, non-contiguous) index i in the matrix of the geminal coefficients
 */
size_t APIGGeminalCoefficients::matrixIndexMajor(size_t vector_index) const {

    return GQCP::matrixIndexMajor(vector_index, this->K);
}


/**
 *  @param vector_index     the vector index of the geminal coefficient
 *
 *  @return the minor (orbital, superscript, contiguous) index p in the matrix of the geminal coefficients
 */
size_t APIGGeminalCoefficients::matrixIndexMinor(size_t vector_index) const {

    return GQCP::matrixIndexMinor(vector_index, this->K);
}


/**
 *  @param i        the major (geminal, subscript, non-contiguous) index
 *  @param p        the minor (orbital, superscript, contiguous) index
 *
 *  @return the vector index of the geminal coefficient G_i^p
 */
size_t APIGGeminalCoefficients::vectorIndex(size_t i, size_t p) const {

    if (i >= N_P) {
        throw std::invalid_argument("The major index i (subscript) must be smaller than N_P.");
    }

    return GQCP::vectorIndex(i, p, this->K);
}


/**
 *  @param onv      the ONV that is being projected on
 *
 *  @return the overlap of the APIG wave function with the given on, i.e. the projection of the APIG wave function onto that ONV
 */
double APIGGeminalCoefficients::overlap(const ONV& onv) const {

    // For a pure APIG, the formula has to be the most general one


    // Construct the matrix G(m) which only has the occupied columns of G in the given ONV m
    Eigen::MatrixXd G = this->asMatrix();  // geminal coefficients as a matrix
    Eigen::MatrixXd Gm = Eigen::MatrixXd::Zero(this->N_P, this->N_P);

    // TODO: wait until the syntax G(Eigen::placeholders::all, occupation_indices) is released in a stable Eigen release
    for (size_t e = 0; e < this->N_P ; e++) {  // loop over all electrons
        size_t occupation_index = onv.get_occupied_index(e);

        Gm.col(e) = G.col(occupation_index);
    }


    // Calculate the permanent of Gm to obtain the coefficient
    return permanent_ryser(Gm);
}


}  // namespace GQCP
