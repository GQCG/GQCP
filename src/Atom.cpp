#include "Atom.hpp"

#include <cmath>
#include <iomanip>

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
 *  OPERATORS
 */

/**
 *  @return if this is equal to @param other, within the @member tolerance_for_comparison for the coordinates
 */
bool Atom::operator==(const GQCG::Atom& other) const {

    return (this->atomic_number == other.atomic_number) &&
           (std::abs(this->x - other.x) < Atom::tolerance_for_comparison) &&
           (std::abs(this->y - other.y) < Atom::tolerance_for_comparison) &&
           (std::abs(this->z - other.z) < Atom::tolerance_for_comparison);
}


/**
 *  @return if this is smaller than @param other, within the @member tolerance_for_comparison for the coordinates
 *
 *  @member atomic_number takes precedence over @member x, over @member y, over @member z
 */
bool Atom::operator<(const GQCG::Atom& other) const {


    if (this->atomic_number < other.atomic_number) {
        return true;
    }

    if (std::abs(this->x - other.x) < Atom::tolerance_for_comparison) {
        return true;
    }

    if (std::abs(this->y - other.y) < Atom::tolerance_for_comparison) {
        return true;
    }

    return (std::abs(this->z - other.z) < Atom::tolerance_for_comparison);  // small simplification of logic expression at the end
}


/**
 *  Overloading of operator<< for a GQCG::Atom to be used with streams
 */
std::ostream& operator<<(std::ostream& os, const GQCG::Atom& atom) {
    os << std::left << std::setw(3) << GQCG::elements::atomicNumberToElement(atom.atomic_number) << '(' << atom.x << ", " << atom.y << ", " << atom.z << ")\n";
    return os;
}



}  // namespace GQCG
