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
#include "RDM/OneRDM.hpp"


namespace GQCP {



/*
 *  CONSTRUCTORS
 */
/**
 *  @param D    the explicit matrix representation of the 1-RDM
 */
OneRDM::OneRDM(const Eigen::MatrixXd& D) :
    BaseRDM(D.cols()),
    D (D)
{
    // Check if the 1-RDM is represented as a square matrix
    if (D.cols() != D.rows()) {
        throw std::invalid_argument("1-RDMs have to be represented as a square matrix.");
    }
}


/*
 *  OPERATORS
 */
/**
 *  @param other    the other OneRDM
 *
 *  @return if the matrix representation of this 1-RDM is equal to the matrix representation of the other, within the default tolerance specified by isEqualTo()
 */
bool OneRDM::operator==(const OneRDM& other) const {
    return this->isEqualTo(other);
}



/*
 *  PUBLIC METHODS
 */

/**
 *  @param other        the other OneRDM
 *  @param tolerance    the tolerance for equality of the matrix representations
 *
 *  @return if the matrix representation of this 1-RDM is equal to the matrix representation of the other, given a tolerance
 */
bool OneRDM::isEqualTo(const OneRDM& other, double tolerance) const {
    return this->D.isApprox(other.D, tolerance);
}


/**
 *  @return the 1-RDM's trace
 */
double OneRDM::trace() const {
    return this->D.trace();
}


/**
 *  In-place diagonalize the 1-RDM
 *
 *  @return the eigenvectors of the 1-RDM
 */
Eigen::MatrixXd OneRDM::diagonalize() {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes (this->D);

    this->D = saes.eigenvalues().asDiagonal();
    return saes.eigenvectors();
}


}  // namespace GQCP
