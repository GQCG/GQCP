// This file is part of GQCG-GQCP.
//
// Copyright (C) 2017-2020  the GQCG developers
//
// GQCG-GQCP is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GQCG-GQCP is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-GQCP.  If not, see <http://www.gnu.org/licenses/>.

#pragma once


#include "Mathematical/Representation/Matrix.hpp"
#include "Utilities/aliases.hpp"
#include "Utilities/type_traits.hpp"

#include <boost/format.hpp>


namespace GQCP {


/*
 *  Forward declaration of ScalarFunction to be used in ScalarFunctionProduct.
 */
template <typename _Valued, typename _Scalar, int _Cols>
class ScalarFunction;


/*
 *  Implementation of ScalarFunctionProduct.
 */

/**
 *  A class template that represents a product of scalar functions (of the same type), such that the evaluation of a product is the product of the evaluations
 *
 *  @tparam T1      the left-hand side scalar function type
 *  @tparam T2      the right-hand side scalar function type
 */
template <typename T1, typename T2 = T1>
class ScalarFunctionProduct:
    public ScalarFunction<product_t<typename T1::Valued, typename T2::Valued>, typename T1::Scalar, T1::Cols> {

    static_assert(std::is_base_of<ScalarFunction<typename T1::Valued, typename T1::Scalar, T1::Cols>, T1>::value, "ScalarFunctionProduct: T1 must inherit from ScalarFunction");
    static_assert(std::is_base_of<ScalarFunction<typename T2::Valued, typename T2::Scalar, T2::Cols>, T2>::value, "ScalarFunctionProduct: T2 must inherit from ScalarFunction");

    static_assert(std::is_same<typename T1::Scalar, typename T2::Scalar>::value, "ScalarFunctionProduct: T1 and T2 should have the same Scalar type");
    static_assert(T1::Cols == T2::Cols, "ScalarFunctionProduct: T1 and T2 should have the same Cols");


private:
    T1 m_lhs;  // left-hand side
    T2 m_rhs;  // right-hand side


public:
    using Valued = product_t<typename T1::Valued, typename T2::Valued>;
    using Scalar = typename T1::Scalar;  // equal to T2::Scalar

    static const auto Cols = T1::Cols;  // equal to T2::Cols


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param lhs      the left-hand side of the product
     *  @param rhs      the right-hand side of the product
     */
    ScalarFunctionProduct(const T1& lhs, const T2& rhs) :
        m_lhs {lhs},
        m_rhs {rhs} {}


    /**
     *  The default constructor.
     */
    ScalarFunctionProduct() = default;


    /**
     *  The constructor for a 'zero' instance given the '0' integer literal.
     */
    ScalarFunctionProduct(const int literal) :
        ScalarFunctionProduct() {

        if (literal != 0) {
            throw std::invalid_argument("ScalarFunctionProduct(const int): Can't convert a non-zero integer to a 'zero' instance.");
        }
    }


    /*
     *  OPERATORS
     */


    /**
     *  @return a textual description of self
     */
    std::string description() const {
        return (boost::format("((%s) * (%s))") % this->m_lhs.description() % this->m_rhs.description()).str();
    }

    /**
     *  @return the left-hand side of this product
     */
    const T1& lhs() const { return this->m_lhs; }

    /**
     *  @param x        the vector/point at which the scalar function product is to be evaluated
     *
     *  @return the product of the evaluated left-hand and right-hand side scalar functions
     */
    Valued operator()(const Vector<Scalar, Cols>& x) const override {
        return this->m_lhs(x) * this->m_rhs(x);
    }


    /**
     *  @return the right-hand side of this product
     */
    const T2& rhs() const { return this->m_rhs; }
};


/*
 *  Implementation of ScalarFunction
 */

/**
 *  A class template representing a mathematical scalar function through overriding operator()
 *
 *  @tparam _Valued     the return type of the scalar function
 *  @tparam _Scalar     the type of the scalars of the input vector
 *  @tparam _Cols       the dimension of the input vector: an integer, or Eigen::Dynamic representing an unknown number of columns at compile time
 */
template <typename _Valued, typename _Scalar, int _Cols>
class ScalarFunction {
public:
    using Valued = _Valued;
    using Scalar = _Scalar;

    static const auto Cols = _Cols;


public:
    /*
     *  DESTRUCTOR
     */
    virtual ~ScalarFunction() = default;


    /*
     *  OPERATORS
     */

    /**
     *  @param x        the vector/point at which the scalar function is to be evaluated
     *
     *  @return the scalar function value at the given point
     */
    virtual _Valued operator()(const Vector<_Scalar, _Cols>& x) const = 0;
};


/*
 *  Aliases related to ScalarFunction.
 */

/**
 *  A SFINAE expression that checks if the type T is a scalar function, i.e. if it derives from ScalarFunction.
 */
template <typename T>
using IsScalarFunction = enable_if_t<std::is_base_of<ScalarFunction<typename T::Valued, typename T::Scalar, T::Cols>, T>::value>;


/**
 *  Multiply one scalar function by another.
 * 
 *  @tparam SF1             the type of the first scalar function
 *  @tparam SF2             the type of the second scalar function
 * 
 *  @note This function is only enabled for actual scalar functions, i.e. functions that derive from ScalarFunction.
 */
template <typename SF1, typename SF2, typename = IsScalarFunction<SF1>, typename = IsScalarFunction<SF2>>
ScalarFunctionProduct<SF1, SF2> operator*(const SF1& lhs, const SF2& rhs) {

    return ScalarFunctionProduct<SF1, SF2>(lhs, rhs);
}


}  // namespace GQCP
