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
#include "HoppingMatrix.hpp"

#include "utilities/miscellaneous.hpp"
#include <iostream>


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param H        the Hubbard hopping matrix
 */
HoppingMatrix::HoppingMatrix(const Eigen::MatrixXd& H) :
    K (H.cols()),
    H (H)
{
    if (H.cols() != H.rows()) {
        throw std::invalid_argument("The given hopping matrix must be square. Did you maybe pass an upper triangle instead of a full hopping matrix?");
    }

    if (!H.transpose().isApprox(H)) {
        throw std::invalid_argument("The given hopping matrix must be symmetric.");
    }
}


/**
 *  Generate the Hubbard hopping matrix from an adjacency matrix and parameters U and t
 *
 *  @param A        the Hubbard adjacency matrix, specifying the connectivity of the Hubbard lattice
 *  @param t        the Hubbard parameter t. Note that a positive value for t means a negative neighbour hopping term
 *  @param U        the Hubbard parameter U
 */
HoppingMatrix::HoppingMatrix(const Eigen::MatrixXd& A, double t, double U) :
    K (A.cols()),
    H (U * Eigen::MatrixXd::Identity(this->K, this->K) - t * A)
{
    if (A.cols() != A.rows()) {
        throw std::invalid_argument("The given adjacency matrix must be square.");
    }

    if (!A.transpose().isApprox(A)) {
        throw std::invalid_argument("The given adjacency matrix must be symmetric.");
    }
}



/*
 *  NAMED CONSTRUCTORS
 */

/**
 *  @param upper_triangle       the upper triangle (in column-major ordering) that specifies the Hubbard hopping matrix
 *
 *  @return the hopping matrix that corresponds to the given upper triangle
 */
HoppingMatrix HoppingMatrix::FromUpperTriangle(const Eigen::VectorXd& upper_triangle) {

    // Check if the given upper triangle is valid
    size_t x = upper_triangle.size();
    size_t K = (static_cast<size_t>(std::sqrt(1 + 8*x) - 1))/2;  // number of rows and columns

    if (K * (K+1) != 2*x) {
        throw std::invalid_argument("The given upper triangle is not a valid upper triangle");
    }


    return HoppingMatrix(fromUpperTriangle(upper_triangle));
}


/**
 *  @param K        the number of lattice sites
 *
 *  @return a random hopping matrix with elements distributed uniformly in [-1.0, 1.0]
 */
HoppingMatrix HoppingMatrix::Random(size_t K) {

    Eigen::VectorXd v = Eigen::VectorXd::Random(K*(K+1)/2);  // random free variables

    return HoppingMatrix::FromUpperTriangle(v);
}



}  // namespace GQCP
