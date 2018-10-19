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
#include "AP1roG/AP1roGGeminalCoefficients.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */
/**
 *  Default constructor setting everything to zero
 */
AP1roGGeminalCoefficients::AP1roGGeminalCoefficients() :
    N_P (0),
    K (0),
    g (Eigen::VectorXd::Zero(0))
{}


/**
 *  Constructor based on given geminal coefficients @param g, which are in vector (row-major) form
 */
AP1roGGeminalCoefficients::AP1roGGeminalCoefficients(const Eigen::VectorXd& g, size_t N_P, size_t K) :
    N_P (N_P),
    K (K),
    g (g)
{
    if (AP1roGGeminalCoefficients::numberOfGeminalCoefficients(N_P, K) != g.size()) {
        throw std::invalid_argument("The specified N_P and K are not compatible with the given vector of geminal coefficients.");
    }
}


/**
 *  Constructor setting the geminal coefficients to zero, based on the number of orbitals @param K and number of electron pairs @param N_P
 */
AP1roGGeminalCoefficients::AP1roGGeminalCoefficients(size_t N_P, size_t K) :
    AP1roGGeminalCoefficients(Eigen::VectorXd::Zero(AP1roGGeminalCoefficients::numberOfGeminalCoefficients(N_P, K)), N_P, K)
{}



/*
 *  PUBLIC METHODS
 */
/**
 *  Construct and @return the geminal coefficients in matrix form
 */
Eigen::MatrixXd AP1roGGeminalCoefficients::asMatrix() const {

    // Initialize the geminal coefficient matrix
    Eigen::MatrixXd G = Eigen::MatrixXd::Zero(this->N_P, this->K);

    // The AP1roG coefficients are the identity matrix in the leftmost (N_P x N_P)-block
    G.topLeftCorner(this->N_P, this->N_P) = Eigen::MatrixXd::Identity(this->N_P, this->N_P);

    // Set the right AP1roG coefficient block
    using RowMajorMatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    Eigen::RowVectorXd g_row = this->g;
    RowMajorMatrixXd B = Eigen::Map<RowMajorMatrixXd, Eigen::RowMajor>(g_row.data(), this->N_P, this->K-this->N_P);
    G.topRightCorner(this->N_P, this->K-this->N_P) = B;

    return G;
}


/*
 *  OPERATORS
 */
/**
 *  @return the geminal coefficient g_mu, in which @param mu is a 'vector index'
 */
double AP1roGGeminalCoefficients::operator()(size_t mu) const {
    return this->g(mu);
}

/**
 *  @ return the geminal coefficient G_i^a, in which @param i is the major index (changes in i are not contiguous) and @param a is the minor index (changes in a are contiguous)
 */
double AP1roGGeminalCoefficients::operator()(size_t i, size_t a) const {
    size_t mu = this->vectorIndex(i, a);
    return this->operator()(mu);
}


/**
 *  @return the number of free geminal coefficients given a number of geminals @param N_P and a number of spatial orbitals @param K
 */
size_t AP1roGGeminalCoefficients::numberOfGeminalCoefficients(size_t N_P, size_t K) {

    // Check if we can have N_P geminals in K orbitals
    if (N_P >= K) {
        throw std::invalid_argument("Can't have that many geminals in this few number of orbitals.");
    }

    return N_P * (K - N_P);
}


/**
 *  For a geminal coefficient g_mu, return its major index in the matrix of geminal coefficients.
 *
 *      Note that:
 *          - the major index is i (i.e. the subscript), since changes in i are not contiguous
 *          - i is in [0 ... N_P[
 */
size_t AP1roGGeminalCoefficients::matrixIndexMajor(size_t vector_index) const {

    return vector_index / (this->K - this->N_P);  // in the mathematical notes, we use the floor function, which is the same as integer division
}


/**
 *  For a geminal coefficient g_mu, return its minor index in the matrix of geminal coefficients.
 *
 *      Note that:
 *          - the minor index is a (i.e. the superscript), since changes in a are contiguous
 *          - a is in [N_P ... K[
 */
size_t AP1roGGeminalCoefficients::matrixIndexMinor(size_t vector_index) const {

    return vector_index % (this->K - this->N_P) + this->N_P;  // we add N_P since we want a to be in [N_P ... K[
}


/**
 *  For a geminal coefficient G_i^a, return its index in the vector of geminal coefficients.
 *
 *      Note that
 *          - i is in [0 ... N_P[       is the 'major' index (i.e. changes in i are not contiguous)
 *          - a is in [N_P ... K[       is the 'minor' index (i.e. changes in a are contiguous)
 */
size_t AP1roGGeminalCoefficients::vectorIndex(size_t i, size_t a) const {

    // Check for invalid values for i and a
    if (i >= this->N_P) {
        throw std::invalid_argument("The major index i (subscript) must be smaller than N_P.");
    }
    if (a < this->N_P) {
        throw std::invalid_argument("The minor index a (superscript) must be larger than or equal to N_P.");
    }


    // The conversion from i and a to a single vector index is just a little more complicated than row-major storage.
    // If we were to use the row-major storage formula, we would end up with
    //      mu = a + (K - N_P) * i
    // but since we would really like our indices abcd (virtual orbitals) to start at N_P, we should subtract N_P accordingly
    return (a - this->N_P) + (this->K - this->N_P) * i;
}


}  // namespace GQCP
