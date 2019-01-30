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
#include "geminals/AP1roGVariables.hpp"

#include "utilities/miscellaneous.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  Default constructor setting everything to zero
 */
AP1roGVariables::AP1roGVariables() :
    BaseAPIGVariables()
{}

    /**
     *  @param x        the variables in a vector representation that is in row-major storage
     *
     *  @param N_P      the number of electron pairs (= the number of geminals)
     *  @param K        the number of spatial orbitals
     */
AP1roGVariables::AP1roGVariables(const Eigen::VectorXd& x, size_t N_P, size_t K) :
    BaseAPIGVariables(x, N_P, K)
{
    if (AP1roGVariables::numberOfVariables(N_P, K) != x.size()) {
        throw std::invalid_argument("The specified N_P and K are not compatible with the given vector of variables.");
    }
}


/**
 *  Constructor that sets the geminal coefficients to zero
 *
 *  @param N_P      the number of electron pairs (= the number of geminals)
 *  @param K        the number of spatial orbitals
 */
AP1roGVariables::AP1roGVariables(size_t N_P, size_t K) :
    AP1roGVariables(Eigen::VectorXd::Zero(AP1roGVariables::numberOfVariables(N_P, K)), N_P, K)
{}



/*
 *  STATIC PUBLIC METHODS
 */

/**
 *  @param N_P      the number of electron pairs (= the number of geminals)
 *  @param K        the number of spatial orbitals
 *
 *  @return the number of 'free' geminal coefficients
 */
size_t AP1roGVariables::numberOfVariables(size_t N_P, size_t K) {

    // Check if we can have N_P geminals in K orbitals
    if (N_P >= K) {
        throw std::invalid_argument("Can't have that many geminals in this few number of orbitals.");
    }

    return N_P * (K - N_P);
}



/*
 *  PUBLIC METHODS
 */

/**
 *  @return the variables in matrix form
 */
Eigen::MatrixXd AP1roGVariables::asMatrix() const {

    using RowMajorMatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    Eigen::RowVectorXd x_row = this->x;
    RowMajorMatrixXd X = Eigen::Map<RowMajorMatrixXd, Eigen::RowMajor>(x_row.data(), this->N_P, this->K-this->N_P);

    return X;
}


/**
 *  @param vector_index     the vector index of the variable
 *
 *  @return the major (geminal, non-contiguous) index i (i.e. the subscript) in the matrix of the variables. Note that i is in [0 ... N_P[
 */
size_t AP1roGVariables::matrixIndexMajor(size_t vector_index) const {

    return GQCP::matrixIndexMajor(vector_index, this->K, this->N_P);
}


/**
 *  @param vector_index     the vector index of the variable
 *
 *  @return the minor (virtual orbital, contiguous) index a (i.e. the subscript) in the matrix of the variables. Note that a is in [N_P ... K[
 */
size_t AP1roGVariables::matrixIndexMinor(size_t vector_index) const {

    return GQCP::matrixIndexMinor(vector_index, this->K, this->N_P);
}


/**
 *  @param i        the major (geminal) index (changes in i are not contiguous)
 *  @param a        the minor (virtual orbital) index (changes in a are contiguous)
 *
 *  @return the vector index of the geminal coefficient G_i^a
 */
size_t AP1roGVariables::vectorIndex(size_t i, size_t a) const {

    if (i >= N_P) {
        throw std::invalid_argument("The major index i (subscript) must be smaller than N_P.");
    }
    if (a < N_P) {
        throw std::invalid_argument("The minor index a (superscript) must be larger than or equal to N_P.");
    }


    return GQCP::vectorIndex(i, a, this->K, this->N_P);
}


}  // namespace GQCP
