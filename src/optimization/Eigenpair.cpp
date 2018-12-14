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
#include "optimization/Eigenpair.hpp"

#include "utilities/linalg.hpp"



namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  A constructor that sets the eigenvalue to zero and the corresponding eigenvector to zeros
 *
 *  @param dimension        the dimension of the eigenvector
 */
Eigenpair::Eigenpair(size_t dimension) :
    eigenvalue (0.0),
    eigenvector (Eigen::VectorXd::Zero(dimension))
{}


/**
 *  @param eigenvalue       the eigenvalue
 *  @param eigenvector      the eigenvector
 */
Eigenpair::Eigenpair(double eigenvalue, const Eigen::VectorXd& eigenvector) :
    eigenvalue (eigenvalue),
    eigenvector (eigenvector)
{}



/*
 *  PUBLIC METHODS
 */

/**
 *  @param other            the other Eigenpair
 *  @param tolerance        a tolerance for comparison
 *
 *  @return if this Eigenpair is equal to the other: if the eigenvalues and eigenvectors are equal given the tolerance
 */
bool Eigenpair::isEqual(const Eigenpair& other, double tolerance) const {


    if (this->eigenvector.size() != other.get_eigenvector().size()) {
        throw std::invalid_argument("Can't compare eigenpairs with eigenvectors of different dimension.");
    }

    if (std::abs(this->eigenvalue - other.eigenvalue) < tolerance) {
        if (GQCP::areEqualEigenvectors(this->eigenvector, other.get_eigenvector(), tolerance)) {
            return true;
        }
    }

    return false;
}


}  // namespace GQCP
