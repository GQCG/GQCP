#ifndef CartesianExponents_hpp
#define CartesianExponents_hpp

#include "CartesianDirection.hpp"

#include <array>


namespace GQCP {


/**
 *  A class that represents Cartesian exponents
 *
 *  Note that we are not deriving from std::array<size_t, 3> because we do not want to support any std::sort algorithms on CartesianExponents: Cartesian components are ordered
 */
class CartesianExponents {
private:
    std::array<size_t, 3> exponents;  // ordered in x, y, z


public:
    // CONSTRUCTORS
    /**
     *  @param exponents        the Cartesian exponents
     */
    CartesianExponents(const std::array<size_t, 3>& exponents);


    // OPERATORS
    /**
     *  @param rhs      the right-hand side of the operator <
     *
     *  @return if these Cartesian exponents are 'smaller' than the ones on the right-hand side. The following logic is used: lhs < rhs
     *      - if lhs's angular momentum is smaller
     *      - if both angular momenta are equal, x takes precedence over y, over z
     *
     *  This means that {1, 0, 0}(x) < {2, 0, 0}(x^2), and {2, 0, 0}(x^2) < {1, 1, 0}(xy)
     */
    bool operator<(const CartesianExponents& rhs) const;

    /**
     *  @param rhs      the right-hand side of the operator ==
     *
     *  @return if the Cartesian exponents are considered equal
     */
    bool operator==(const CartesianExponents& rhs) const;


    // PUBLIC METHODS
    /**
     *  @return the underlying values of the Cartesian exponents
     */
    const std::array<size_t, 3>& values() const;

    /**
     *  @param direction        the direction (x,y,z) whose exponent should be returned
     *
     *  @return the exponent in the given direction
     */
    size_t value(CartesianDirection direction) const;

    /**
     *  @return the angular momentum corresponding to these exponents
     */
    size_t angularMomentum() const;

    /**
     *  @return the exponent belonging to x
     */
    size_t x() const;

    /**
     *  @return a reference to the exponent belonging to x
     */
    size_t& x();

    /**
     *  @return the exponent belonging to y
     */
    size_t y() const;

    /**
     *  @return a reference to the exponent belonging to y
     */
    size_t& y();

    /**
     *  @return the exponent belonging to z
     */
    size_t z() const;

    /**
     *  @return a reference to the exponent belonging to z
     */
    size_t& z();
};


}  // namespace GQCP


#endif  /* CartesianExponents_hpp */
