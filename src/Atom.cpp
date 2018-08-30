#include "Atom.hpp"

#include <cmath>

#include "elements.hpp"


namespace GQCG {


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
 *  PUBLIC METHODS
 */

/**
 *  @return if this is equal to @param other, within the given @param tolerance for the coordinates
 */
bool Atom::isEqualTo(const GQCG::Atom& other, double tolerance) const {
    return (this->atomic_number == other.atomic_number) &&
           (std::abs(this->x - other.x) < tolerance) &&
           (std::abs(this->y - other.y) < tolerance) &&
           (std::abs(this->z - other.z) < tolerance);
}


}  // namespace GQCG


/**
 *  operator<< for a GQCG::Atom
 */
std::ostream& operator<<(std::ostream& os, const GQCG::Atom& atom) {
    os << GQCG::elements::atomicNumberToElement(atom.atomic_number) << ": (" << atom.x << ", " << atom.y << ", " << atom.z << ")\n";
    return os;
}
