// This file is part of GQCG-GQCP.
//
// Copyright (C) 2017-2020  the GQCG developers
//
// GQCG-GQCP is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GQCG-GQCP is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-GQCP.  If not, see <http://www.gnu.org/licenses/>.

#include "Mathematical/Optimization/Eigenproblem/Eigenpair.hpp"

#include "Mathematical/Representation/Matrix.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param eigenvalue       the eigenvalue
 *  @param eigenvector      the eigenvector
 */
Eigenpair::Eigenpair(const double eigenvalue, const VectorX<double>& eigenvector) :
    m_eigenvalue {eigenvalue},
    m_eigenvector {eigenvector} {}


/**
 *  A constructor that sets the eigenvalue to zero and the corresponding eigenvector to zeros
 *
 *  @param dimension        the dimension of the eigenvector
 */
Eigenpair::Eigenpair(const size_t dimension) :
    Eigenpair(0.0, VectorX<double>::Zero(dimension)) {}


/*
 *  PUBLIC METHODS
 */

/**
 *  @param other            the other Eigenpair
 *  @param tolerance        a tolerance for comparison
 *
 *  @return if this Eigenpair is equal to the other: if the eigenvalues and eigenvectors are equal given the tolerance
 */
bool Eigenpair::isEqualTo(const Eigenpair& other, const double tolerance) const {


    if (this->eigenvector().size() != other.eigenvector().size()) {
        throw std::invalid_argument("Eigenpair::isEqualTo(Eigenpair, double): Can't compare eigenpairs with eigenvectors of different dimension.");
    }

    if (std::abs(this->eigenvalue() - other.eigenvalue()) < tolerance) {
        if ((this->eigenvector()).isEqualEigenvectorAs(other.eigenvector(), tolerance)) {
            return true;
        }
    }

    return false;
}


}  // namespace GQCP
