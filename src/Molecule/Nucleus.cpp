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

#include <cmath>
#include <iomanip>

#include "elements.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param atomic_number        the atomic number (Z) of the nucleus
 *  @param x                    the x-position of the nucleus in bohr
 *  @param y                    the y-position of the nucleus in bohr
 *  @param z                    the z-position of the nucleus in bohr
 */
Nucleus::Nucleus(size_t atomic_number, double x, double y, double z) :
    atomic_number (atomic_number)
{
    position << x, y, z;
}


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
 *  @param other        the other nucleus
 *
 *  @return if this nucleus is equal to the other, within a default tolerance for the coordinates
 */
bool Nucleus::operator==(const Nucleus& other) const {

    return this->isEqualTo(other, Nucleus::tolerance_for_comparison);
}


/**
 *  @param other        the other nucleus
 *
 *  @return if this nucleus is not equal to the other, within a default tolerance for the coordinates
 */
bool Nucleus::operator!=(const Nucleus& other) const {
    return !this->operator==(other);
}


/**
 *  A custom implementation for the comparison (and thus ordening) of nuclei. The atomic_number takes precedence over the x-coordinate, which takes precedence over the y-coordinate, which in turn takes precedence over the z-coordinate
 *
 *  @param other        the other nucleus
 *
 *  @return if this nucleus is 'smaller' than the other, within a default tolerance for the coordinates
 */
bool Nucleus::operator<(const Nucleus& other) const {

    return this->isSmallerThan(other, Nucleus::tolerance_for_comparison);
}


/**
 *  Overloading of operator<< for a Nucleus to be used with ostreams
 *
 *  @param os           the output stream to which the nucleus should be concatenated
 *  @param nucleus      the nucleus which should be concatenated to the output stream
 *
 *  @return the updated output stream
 */
std::ostream& operator<<(std::ostream& os, const Nucleus& nucleus) {
    os << std::left << std::setw(3) << elements::atomicNumberToElement(nucleus.atomic_number) << '(' << nucleus.position.x() << ", " << nucleus.position.y() << ", " << nucleus.position.z() << ")\n";
    return os;
}


/*
 *  PUBLIC METHODS
 */

/**
 *  @param other        the other nucleus
 *  @param tolerance    the tolerance for equality of positions
 *
 *  @return if this nucleus is equal to the other
 */
bool Nucleus::isEqualTo(const Nucleus& other, double tolerance) const {

    return (this->atomic_number == other.atomic_number) &&
           (std::abs(this->position.x() - other.position.x()) < tolerance) &&
           (std::abs(this->position.y() - other.position.y()) < tolerance) &&
           (std::abs(this->position.z() - other.position.z()) < tolerance);
}


/**
 *  A custom implementation for the comparison (and thus ordening) of nuclei. The atomic_number takes precedence over the x-coordinate, which takes precedence over the y-coordinate, which in turn takes precedence over the z-coordinate
 *
 *  @param other        the other nucleus
 *  @param tolerance    the tolerance for equality of positions
 *
 *  @return if this nucleus is 'smaller' than the other, within a default tolerance for the coordinates
 */
bool Nucleus::isSmallerThan(const Nucleus& other, double tolerance) const {

    if (this->atomic_number < other.atomic_number) {
        return true;
    } else if (this->atomic_number > other.atomic_number) {
        return false;
    } else {  // the atomic numbers are equal

        if (std::abs(this->position.x() - other.position.x()) > tolerance) {  // the difference is meaningful
            return (this->position.x() < other.position.x());
        } else {  // the x-coordinates are considered equal

            if (std::abs(this->position.y() - other.position.y()) > tolerance) {  // the difference is meaningful
                return (this->position.y() < other.position.y());
            } else {  // the y-coordinates are considered equal

                if (std::abs(this->position.z() - other.position.z()) > tolerance) {  // the difference is meaningful
                    return (this->position.z() < other.position.z());
                } else {  // the z-coordinates are considered equal
                    return false;
                }
            }  // else y
        }  // else x
    }  // else atomic_number
}


/**
 *  @param other        the other nucleus
 *
 *  @return the Euclidian distance between this nucleus and the other
 */
double Nucleus::calculateDistance(const Nucleus& other) const {

    return (this->position - other.position).norm();
}


}  // namespace GQCP
