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
#include "Atom.hpp"

#include <cmath>
#include <iomanip>

#include "elements.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param atomic_number        the atomic number (Z) of the atom
 *  @param x                    the x-position of the atom
 *  @param y                    the y-position of the atom
 *  @param z                    the z-position of the atom
 */
Atom::Atom(size_t atomic_number, double x, double y, double z) :
    atomic_number (atomic_number),
    position (Eigen::Vector3d {x, y, z})
{}



/*
 *  OPERATORS
 */

/**
 *  @param other        the other atom
 *
 *  @return if this atom is equal to the other, within a default tolerance for the coordinates
 */
bool Atom::operator==(const GQCP::Atom& other) const {

    return this->isEqualTo(other, Atom::tolerance_for_comparison);
}


/**
 *  A custom implementation for the comparison (and thus ordening) of atoms. The atomic_number takes precedence over the x-coordinate, which takes precedence over the y-coordinate, which in turn takes precedence over the z-coordinate
 *
 *  @param other        the other atom
 *
 *  @return if this atom is 'smaller' than the other, within a default tolerance for the coordinates
 */
bool Atom::operator<(const GQCP::Atom& other) const {

    return this->isSmallerThan(other, Atom::tolerance_for_comparison);
}


/**
 *  Overloading of operator<< for a GQCP::Atom to be used with ostreams
 *
 *  @param os       the output stream to which the atom should be concatenated
 *  @param atom     the atom which should be concatenated to the output stream
 *
 *  @return the updated output stream
 */
std::ostream& operator<<(std::ostream& os, const GQCP::Atom& atom) {
    os << std::left << std::setw(3) << GQCP::elements::atomicNumberToElement(atom.atomic_number) << '(' << atom.position.x() << ", " << atom.position.y() << ", " << atom.position.z() << ")\n";
    return os;
}


/*
 *  PUBLIC METHODS
 */

/**
 *  @param other        the other atom
 *  @param tolerance    the tolerance for equality of positions
 *
 *  @return if this atom is equal to the other
 */
bool Atom::isEqualTo(const GQCP::Atom& other, double tolerance) const {

    return (this->atomic_number == other.atomic_number) &&
           (std::abs(this->position.x() - other.position.x()) < tolerance) &&
           (std::abs(this->position.y() - other.position.y()) < tolerance) &&
           (std::abs(this->position.z() - other.position.z()) < tolerance);
}


/**
 *  A custom implementation for the comparison (and thus ordening) of atoms. The atomic_number takes precedence over the x-coordinate, which takes precedence over the y-coordinate, which in turn takes precedence over the z-coordinate
 *
 *  @param other        the other atom
 *  @param tolerance    the tolerance for equality of positions
 *
 *  @return if this atom is 'smaller' than the other, within a default tolerance for the coordinates
 */
bool Atom::isSmallerThan(const GQCP::Atom& other, double tolerance) const {

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
 *  @param other        the other atom
 *
 *  @return the Euclidian distance between this atom and the other
 */
double Atom::calculateDistance(const GQCP::Atom& other) const {

    return (this->position - other.position).norm();
}


}  // namespace GQCP
