// This file is part of GQCG-numopt.
// 
// Copyright (C) 2017-2018  the GQCG developers
// 
// GQCG-numopt is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// GQCG-numopt is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-numopt.  If not, see <http://www.gnu.org/licenses/>.
#include "Eigenpair.hpp"

#include <cpputil.hpp>



namespace numopt {
namespace eigenproblem {


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
bool Eigenpair::isEqual(const numopt::eigenproblem::Eigenpair& other, double tolerance) const {


    if (this->eigenvector.size() != other.get_eigenvector().size()) {
        throw std::invalid_argument("Can't compare eigenpairs with eigenvectors of different dimension.");
    }

    if (std::abs(this->eigenvalue - other.eigenvalue) < tolerance) {
        if (cpputil::linalg::areEqualEigenvectors(this->eigenvector, other.get_eigenvector(), tolerance)) {
            return true;
        }
    }

    return false;
}


}  // namespace eigenproblem
}  // namespace numopt
