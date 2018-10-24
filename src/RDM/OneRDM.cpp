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

OneRDM::OneRDM(const Eigen::MatrixXd& D) :
    BaseRDM (D.cols()),
    D (D)
{
    // Check if the 1-RDM is represented as a square matrix
    if (D.cols() != D.rows()) {
        throw std::invalid_argument("1-RDMs have to be represented as a square matrix.");
    }
}



/*
 *  PUBLIC METHODS
 */

/**
 *  @return the 1-RDM's trace
 */
double OneRDM::trace(){
    return this->D.trace();
}


/**
 *  diagonalizes the 1-RDM and @returns the eigenvectors
 */
Eigen::MatrixXd OneRDM::diagonalize() {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes (this->D);
    this->D = saes.eigenvalues().asDiagonal();
    return saes.eigenvectors();
}


}  // namespace GQCP
