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
#include "Molecule/Nucleus.hpp"

#include "elements.hpp"

#include <cmath>
#include <iomanip>


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param Z        the atomic number (Z) of the nucleus
 *  @param x        the x-position of the nucleus in bohr
 *  @param y        the y-position of the nucleus in bohr
 *  @param z        the z-position of the nucleus in bohr
 */
Nucleus::Nucleus(const size_t Z, const double x, const double y, const double z) :
    Z (Z)
{
    this->R.setZero();
    this->R << x, y, z;
}


/**
 *  @param Z            the atomic number (Z) of the nucleus
 *  @param position     the position of the nucleus in bohr
 */
Nucleus::Nucleus(const size_t Z, const Vector<double, 3>& position) :
    Nucleus(Z, position.x(), position.y(), position.z())
{}


/**
 *  Default constructor, creating a 'ghost' nucleus (i.e. Bq) in the origin
 */
Nucleus::Nucleus() :
    Nucleus(0,  0.0, 0.0, 0.0)  // Z = 0
{}



/*
 *  OPERATORS
 */

/**
 *  Overloading of operator<< for a Nucleus to be used with ostreams
 *
 *  @param os           the output stream to which the nucleus should be concatenated
 *  @param nucleus      the nucleus which should be concatenated to the output stream
 *
 *  @return the updated output stream
 */
std::ostream& operator<<(std::ostream& os, const Nucleus& nucleus) {
    os << std::left << std::setw(3) << elements::atomicNumberToElement(nucleus.charge()) << '(' << nucleus.position().x() << ", " << nucleus.position().y() << ", " << nucleus.position().z() << ")\n";
    return os;
}



/*
 *  STATIC PUBLIC METHODS
 */

/**
 *  @return a functor that can be used in sorting atoms. It features a custom implementation, in which the x-coordinate takes precedence over the y-coordinate, which in turn takes precedence over the z-coordinate
 */
std::function<bool(const Nucleus&, const Nucleus&)> Nucleus::sortComparer(const double tolerance) {

    const auto sort_comparer_lambda = [tolerance] (const Nucleus& lhs, const Nucleus& rhs) {

        if (std::abs(lhs.position().x() - rhs.position().x()) > tolerance) {  // the difference is meaningful
            return (lhs.position().x() < rhs.position().x());
        } else {  // the x-coordinates are considered equal

            if (std::abs(lhs.position().y() - rhs.position().y()) > tolerance) {  // the difference is meaningful
                return (lhs.position().y() < rhs.position().y());
            } else {  // the y-coordinates are considered equal

                if (std::abs(lhs.position().z() - rhs.position().z()) > tolerance) {  // the difference is meaningful
                    return (lhs.position().z() < rhs.position().z());
                } else {  // the z-coordinates are considered equal
                    return false;
                }
            }  // else y
        }  // else x
    };

    return sort_comparer_lambda;
}


/**
 *  @return a functor that can be used in checking atom for equality. Atoms are equal if their charge and position are equal
 */
std::function<bool(const Nucleus&, const Nucleus&)> Nucleus::equalityComparer(const double tolerance) {

    const auto equality_comparer_lambda = [tolerance] (const Nucleus& lhs, const Nucleus& rhs) {

        return (lhs.charge() == rhs.charge()) &&
           (std::abs(lhs.position().x() - rhs.position().x()) < tolerance) &&
           (std::abs(lhs.position().y() - rhs.position().y()) < tolerance) &&
           (std::abs(lhs.position().z() - rhs.position().z()) < tolerance);
    };

    return equality_comparer_lambda;
}


/*
 *  PUBLIC METHODS
 */

/**
 *  @param other        the other nucleus
 *
 *  @return the Euclidian distance between this nucleus and the other
 */
double Nucleus::calculateDistance(const Nucleus& other) const {

    return (this->position() - other.position()).norm();
}


}  // namespace GQCP
