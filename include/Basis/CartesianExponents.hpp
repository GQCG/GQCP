#ifndef CartesianExponents_hpp
#define CartesianExponents_hpp

#include "CartesianDirection.hpp"

#include <array>


namespace GQCP {




//    typedef array __self;
//    typedef _Tp                                   value_type;
//    typedef value_type&                           reference;
//    typedef const value_type&                     const_reference;
//    typedef value_type*                           iterator;
//    typedef const value_type*                     const_iterator;
//    typedef value_type*                           pointer;
//    typedef const value_type*                     const_pointer;
//    typedef size_t                                size_type;
//    typedef ptrdiff_t                             difference_type;
//    typedef std::reverse_iterator<iterator>       reverse_iterator;
//    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
//
//    value_type __elems_[_Size > 0 ? _Size : 1];
//
//    // No explicit construct/copy/destroy for aggregate type
//    _LIBCPP_INLINE_VISIBILITY void fill(const value_type& __u)
//    {_VSTD::fill_n(__elems_, _Size, __u);}
//    _LIBCPP_INLINE_VISIBILITY
//    void swap(array& __a) _NOEXCEPT_(_Size == 0 || __is_nothrow_swappable<_Tp>::value)
//    { __swap_dispatch((std::integral_constant<bool, _Size == 0>()), __a); }
//
//    _LIBCPP_INLINE_VISIBILITY
//    void __swap_dispatch(std::true_type, array&) {}
//
//    _LIBCPP_INLINE_VISIBILITY
//    void __swap_dispatch(std::false_type, array& __a)
//    { _VSTD::swap_ranges(__elems_, __elems_ + _Size, __a.__elems_);}
//
//    // iterators:
//    _LIBCPP_INLINE_VISIBILITY _LIBCPP_CONSTEXPR_AFTER_CXX14
//    iterator begin() _NOEXCEPT {return iterator(__elems_);}
//    _LIBCPP_INLINE_VISIBILITY _LIBCPP_CONSTEXPR_AFTER_CXX14
//    const_iterator begin() const _NOEXCEPT {return const_iterator(__elems_);}
//    _LIBCPP_INLINE_VISIBILITY _LIBCPP_CONSTEXPR_AFTER_CXX14
//    iterator end() _NOEXCEPT {return iterator(__elems_ + _Size);}
//    _LIBCPP_INLINE_VISIBILITY _LIBCPP_CONSTEXPR_AFTER_CXX14
//    const_iterator end() const _NOEXCEPT {return const_iterator(__elems_ + _Size);}

/**
 *  A class that represents Cartesian exponents
 *
 *  Note that we are not deriving from std::array<size_t, 3> because we want to overwrite the behavior of operator<
 */
class CartesianExponents {
private:
    std::array<size_t, 3> exponents;  // ordered in x, y, z


public:
    using difference_type = std::ptrdiff_t;
    using value_type = size_t;
    using pointer = value_type*;
    using reference = value_type&;
    using iterator_category = std::forward_iterator_tag;

    using iterator = value_type*;


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
    const std::array<value_type, 3>& values() const;

    /**
     *  @param direction        the direction (x,y,z) whose exponent should be returned
     *
     *  @return the exponent in the given direction
     */
    value_type value(CartesianDirection direction) const;

    /**
     *  @return the exponent belonging to x
     */
    value_type x() const;

    /**
     *  @return a reference to the exponent belonging to x
     */
    reference x();

    /**
     *  @return the exponent belonging to y
     */
    value_type y() const;

    /**
     *  @return a reference to the exponent belonging to y
     */
    reference y();

    /**
     *  @return the exponent belonging to z
     */
    value_type z() const;

    /**
     *  @return a reference to the exponent belonging to z
     */
    reference z();

    /**
     *  @return the angular momentum corresponding to these exponents
     */
    value_type angularMomentum() const;

    /**
     *  @return an iterator to this' begin
     */
    iterator begin();

    /**
     *  @return an iterator to this' end
     */
    iterator end();
};


}  // namespace GQCP


#endif  /* CartesianExponents_hpp */
