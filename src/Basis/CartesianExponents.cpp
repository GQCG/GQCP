#include "CartesianExponents.hpp"

#include <numeric>


namespace GQCP {

/*
 *  CONSTRUCTORS
 */

/**
 *  @param exponents        the Cartesian exponents
 */
CartesianExponents::CartesianExponents(const std::array<size_t, 3>& exponents) :
    exponents (exponents)
{}



/*
 *  OPERATORS
 */

/**
 *  @param rhs      the right-hand side of the operator <
 *
 *  @return if these Cartesian exponents are 'smaller' than the ones on the right-hand side. The following logic is used: lhs < rhs
 *      - if lhs's angular momentum is smaller
 *      - if both angular momenta are equal, a larger x takes precedence over a larger y, over a larger z
 */
bool CartesianExponents::operator<(const CartesianExponents& rhs) const {

    // Compare angular momenta
    if (this->angularMomentum() < rhs.angularMomentum()) { return true; }
    else if (this->angularMomentum() > rhs.angularMomentum()) { return false; }

    else {
        // Compare x-exponents
        if (this->x() > rhs.x()) { return true; }
        else if (this->x() < rhs.x()) { return false; }

        else {
            // Compare y-exponents
            if (this->y() > rhs.y()) { return true; }
            else if (this->y() < rhs.y()) { return false; }

            else {
                // Compare z-exponents
                if (this->z() >= rhs.z()) { return true; }  // if, at the end, the z exponents are equal, then the Cartesian exponents should be considered equal
                else { return false; }
            }  // same y
        }  // same x
    }  // same angular momentum
}


/**
 *  @param rhs      the right-hand side of the operator ==
 *
 *  @return if the Cartesian exponents are considered equal
 */
bool CartesianExponents::operator==(const CartesianExponents& rhs) const {
    return this->exponents == rhs.exponents;
}


/*
 *  PUBLIC METHODS
 */

/**
 *  @return the underlying values of the Cartesian components
 */
const std::array<size_t, 3>& CartesianExponents::values() const {
    return this->exponents;
}


/**
 *  @return the angular momentum corresponding to these exponents
 */
size_t CartesianExponents::angularMomentum() const {
    return std::accumulate(this->exponents.begin(), this->exponents.end(), 0);  // 0 is starting value of the sum
}


/**
 *  @return the exponent belonging to x
 */
size_t CartesianExponents::x() const {
    return this->exponents[0];
}


/**
 *  @return a reference to the exponent belonging to x
 */
size_t& CartesianExponents::x() {
    return this->exponents[0];
}


/**
 *  @return the exponent belonging to y
 */
size_t CartesianExponents::y() const {
    return this->exponents[1];
}


/**
 *  @return a reference to the exponent belonging to y
 */
size_t& CartesianExponents::y() {
    return this->exponents[1];
}


/**
 *  @return the exponent belonging to z
 */
size_t CartesianExponents::z() const {
    return this->exponents[2];
}


/**
 *  @return a reference to the exponent belonging to z
 */
size_t& CartesianExponents::z() {
    return this->exponents[2];
}


}  // namespace GQCP
