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
#ifndef GQCP_ONERDM_HPP
#define GQCP_ONERDM_HPP


#include "RDM/BaseRDM.hpp"

#include <Eigen/Dense>


namespace GQCP {

/**
 *  A class that holds the matrix representation of a 1-RDM
 */
class OneRDM : public BaseRDM {
private:
    Eigen::MatrixXd D;


public:
    // CONSTRUCTORS
    explicit OneRDM(const Eigen::MatrixXd& D);


    // OPERATORS
    double operator()(size_t p, size_t q) const { return this->D(p,q); }


    // GETTERS
    Eigen::MatrixXd get_matrix_representation() const { return this->D; }


    // PUBLIC METHODS
    /**
     *  @return the 1-RDM's trace
     */
    double trace();

    /**
     *  diagonalizes the 1-RDM and @returns the eigenvectors
     */
    Eigen::MatrixXd diagonalize();
};


}  // namespace GQCP


#endif  // GQCP_ONERDM_HPP
