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

#pragma once


#include "Mathematical/Representation/Matrix.hpp"


namespace GQCP {


/**
 *  A container class to store an eigenpair, i.e. an eigenvector with a corresponding eigenvalue.
 */
template <typename _Scalar = double>
class Eigenpair {
public:
    using Scalar = _Scalar;

private:
    Scalar m_eigenvalue;
    VectorX<Scalar> m_eigenvector;


public:
    // CONSTRUCTORS

    /**
     *  @param eigenvalue       The eigenvalue.
     *  @param eigenvector      The eigenvector.
     */
    Eigenpair(const Scalar eigenvalue, const VectorX<Scalar>& eigenvector) :
        m_eigenvalue {eigenvalue},
        m_eigenvector {eigenvector} {}

    /**
     *  A constructor that sets the eigenvalue to zero and the corresponding eigenvector to zeros
     *
     *  @param dimension        The dimension of the eigenvector.
     */
    explicit Eigenpair(const size_t dimension = 1) :
        Eigenpair(0.0, VectorX<double>::Zero(dimension)) {}


    // PUBLIC METHODS

    /**
     *  @return The eigenvalue associated to this eigenpair.
     */
    Scalar eigenvalue() const { return this->m_eigenvalue; };

    /**
     * A function that returns the eigenvector associated to this eigenpair.
     * 
     *  @return the eigenvector associated to this eigenpair
     */
    const VectorX<Scalar>& eigenvector() const { return this->m_eigenvector; };

    /**
     *  Check if this Eigenpair is equal to the other Eigenpair.
     * 
     *  @param other            The other Eigenpair.
     *  @param tolerance        A tolerance for comparison.
     *
     *  @return if this Eigenpair is equal to the other: if the eigenvalues and eigenvectors are equal given the tolerance
     */
    bool isEqualTo(const Eigenpair& other, const double tolerance = 1.0e-08) const {

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
};


}  // namespace GQCP
