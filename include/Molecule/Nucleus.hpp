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
#ifndef GQCP_NUCLEUS_HPP
#define GQCP_NUCLEUS_HPP


#include "Mathematical/Matrix.hpp"

#include <functional>
#include <iostream>
#include <stdlib.h>


namespace GQCP {


/**
 *  A class that represents a nucleus: it has a charge and a position in space
 */
class Nucleus {
private:
    size_t Z;  // atomic number, i.e. the charge
    Vector<double, 3> R;  // in bohr


public:
    // CONSTRUCTORS

    /**
     *  @param Z        the atomic number (Z) of the nucleus
     *  @param x        the x-position of the nucleus in bohr
     *  @param y        the y-position of the nucleus in bohr
     *  @param z        the z-position of the nucleus in bohr
     */
    Nucleus(const size_t Z, const double x, const double y, const double z);

    /**
     *  @param Z            the atomic number (Z) of the nucleus
     *  @param position     the position of the nucleus in bohr
     */
    Nucleus(const size_t Z, const Vector<double, 3>& position);

    /**
     *  Default constructor, creating a 'ghost' nucleus (i.e. Bq) in the origin
     */
    Nucleus();


    /*
     *  OPERATORS
     */

    /**
     *  Overloading of operator<< for a Nucleus to be used with ostreams
     *
     *  @param os       the output stream to which the nucleus should be concatenated
     *  @param nucleus     the nucleus which should be concatenated to the output stream
     *
     *  @return the updated output stream
     */
    friend std::ostream& operator<<(std::ostream& os, const Nucleus& nucleus);


    // STATIC PUBLIC METHODS

    /**
     *  @return a functor that can be used in sorting atoms. It features a custom implementation, in which the x-coordinate takes precedence over the y-coordinate, which in turn takes precedence over the z-coordinate
     */
    static std::function<bool(const Nucleus&, const Nucleus&)> sortComparer(const double tolerance=1.0e-08);

    /**
     *  @return a functor that can be used in checking atom for equality. Atoms are equal if their charge and position are equal
     */
    static std::function<bool(const Nucleus&, const Nucleus&)> equalityComparer(const double tolerance=1.0e-08);


    // PUBLIC METHODS

    /**
     *  @return the charge of this nucleus
     */
    size_t charge() const { return this->Z; }

    /**
     *  @return the position of this nucleus
     */
    const Vector<double, 3>& position() const { return this->R; }

    /**
     *  @param other        the other nucleus
     *
     *  @return the Euclidian distance between this nucleus and the other
     */
    double calculateDistance(const Nucleus& other) const;
};


}  // namespace GQCP



#endif  // GQCP_NUCLEUS_HPP
