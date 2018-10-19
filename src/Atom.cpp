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
    x (x),
    y (y),
    z (z)
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
    os << std::left << std::setw(3) << GQCP::elements::atomicNumberToElement(atom.atomic_number) << '(' << atom.x << ", " << atom.y << ", " << atom.z << ")\n";
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
           (std::abs(this->x - other.x) < tolerance) &&
           (std::abs(this->y - other.y) < tolerance) &&
           (std::abs(this->z - other.z) < tolerance);
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

        if (std::abs(this->x - other.x) > tolerance) {  // the difference is meaningful
            return (this->x < other.x);
        } else {  // the x-coordinates are considered equal

            if (std::abs(this->y - other.y) > tolerance) {  // the difference is meaningful
                return (this->y < other.y);
            } else {  // the y-coordinates are considered equal

                if (std::abs(this->z - other.z) > tolerance) {  // the difference is meaningful
                    return (this->z < other.z);
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

    return std::sqrt(std::pow(this->x - other.x, 2)
                     + std::pow(this->y - other.y, 2)
                     + std::pow(this->z - other.z, 2));
}


}  // namespace GQCP
