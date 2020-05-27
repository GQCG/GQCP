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
 *  A container class to store an eigenpair, i.e. an eigenvector with a corresponding eigenvalue
 */
class Eigenpair {
private:
    double m_eigenvalue;
    VectorX<double> m_eigenvector;


public:
    // CONSTRUCTORS

    /**
     *  @param eigenvalue       the eigenvalue
     *  @param eigenvector      the eigenvector
     */
    Eigenpair(const double eigenvalue, const VectorX<double>& eigenvector);

    /**
     *  A constructor that sets the eigenvalue to zero and the corresponding eigenvector to zeros
     *
     *  @param dimension        the dimension of the eigenvector
     */
    Eigenpair(const size_t dimension = 1);


    // PUBLIC METHODS

    /**
     *  @return the eigenvalue associated to this eigenpair
     */
    double eigenvalue() const { return this->m_eigenvalue; };

    /**
     *  @return the eigenvector associated to this eigenpair
     */
    const VectorX<double>& eigenvector() const { return this->m_eigenvector; };

    /**
     *  @param other            the other Eigenpair
     *  @param tolerance        a tolerance for comparison
     *
     *  @return if this Eigenpair is equal to the other: if the eigenvalues and eigenvectors are equal given the tolerance
     */
    bool isEqualTo(const Eigenpair& other, const double tolerance = 1.0e-08) const;
};


}  // namespace GQCP
