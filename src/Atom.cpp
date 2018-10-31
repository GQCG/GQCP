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
 *  Constructor based on a given @param atomic_number and the coordinates @param x, @param y, @param z
 */
Atom::Atom(size_t atomic_number, double x, double y, double z) :
    atomic_number (atomic_number),
    position (Eigen::Vector3d {x, y, z})
{}



/*
 *  OPERATORS
 */

/**
 *  @return if this is equal to @param other, within the @member tolerance_for_comparison for the coordinates
 */
bool Atom::operator==(const GQCP::Atom& other) const {

    return this->isEqualTo(other, Atom::tolerance_for_comparison);
}


/**
 *  @return if this is smaller than @param other, within the @member tolerance_for_comparison for the coordinates
 *
 *  @member atomic_number takes precedence over @member x, over @member y, over @member z
 */
bool Atom::operator<(const GQCP::Atom& other) const {

    return this->isSmallerThan(other, Atom::tolerance_for_comparison);
}


/**
 *  Overloading of operator<< for a GQCP::Atom to be used with streams
 */
std::ostream& operator<<(std::ostream& os, const GQCP::Atom& atom) {
    os << std::left << std::setw(3) << GQCP::elements::atomicNumberToElement(atom.atomic_number) << '(' << atom.position.x() << ", " << atom.position.y() << ", " << atom.position.z() << ")\n";
    return os;
}


/*
 *  PUBLIC METHODS
 */

/**
 *  @return if this is equal to @param other, within the given @param tolerance for the coordinates
 */
bool Atom::isEqualTo(const GQCP::Atom& other, double tolerance) const {

    return (this->atomic_number == other.atomic_number) &&
           (std::abs(this->position.x() - other.position.x()) < tolerance) &&
           (std::abs(this->position.y() - other.position.y()) < tolerance) &&
           (std::abs(this->position.z() - other.position.z()) < tolerance);
}


/**
 *  @return if this is smaller than @param other, within the given @param tolerance for the coordinates
 *
 *  @member atomic_number takes precedence over @member x, over @member y, over @member z
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
 * @return the distance between this and @param other
 */
double Atom::calculateDistance(const GQCP::Atom& other) const {

    return (this->position - other.position).norm();
}


}  // namespace GQCP
