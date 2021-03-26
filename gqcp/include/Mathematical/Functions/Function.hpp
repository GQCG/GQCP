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


#include "Utilities/aliases.hpp"
#include "Utilities/type_traits.hpp"

#include <boost/format.hpp>


namespace GQCP {


/*
 *  MARK: Forward declarations
 */
template <typename OutputType, typename InputType>
class Function;


/*
 *  MARK: `FunctionProduct` implementation
 */

/**
 *  A template that represents a product of functions. In this way, the evaluation of a function product is the product of the respective evaluations.
 *
 *  @tparam T1          The left-hand side function type.
 *  @tparam T2          The right-hand side function type
 */
template <typename T1, typename T2 = T1>
class FunctionProduct:
    public Function<product_t<typename T1::OutputType, typename T2::OutputType>, typename T1::InputType> {

    static_assert(std::is_same<typename T1::InputType, typename T2::InputType>::value, "FunctionProduct: T1 and T2 should have the same InputType.");

private:
    // The left-hand side of the product.
    T1 m_lhs;

    // The right-hand side of the product.
    T2 m_rhs;


public:
    // The combined output type of the `operator()`.
    using OutputType = product_t<typename T1::OutputType, typename T2::OutputType>;

    // The input type of the `operator()`.
    using InputType = typename T1::InputType;  // Due to the `static_assert`, this is equal to `T2::InputType`.


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  @param lhs          The left-hand side of the product.
     *  @param rhs          The right-hand side of the product.
     */
    FunctionProduct(const T1& lhs, const T2& rhs) :
        m_lhs {lhs},
        m_rhs {rhs} {}


    /**
     *  The default constructor.
     */
    FunctionProduct() = default;


    /**
     *  Construct a 'zero function product' given the '0' integer literal.
     * 
     *  @param zero         The '0' integer literal.
     */
    FunctionProduct(const int zero) :
        FunctionProduct() {

        if (zero != 0) {
            throw std::invalid_argument("FunctionProduct(const int): Can't convert a non-zero integer to a 'zero' instance.");
        }
    }


    /*
     *  MARK: General information
     */

    /**
     *  @return A textual description of this function.
     */
    std::string description() const {
        return (boost::format("((%s) * (%s))") % this->m_lhs.description() % this->m_rhs.description()).str();
    }


    /*
     *  MARK: Access
     */

    /**
     *  @return The left-hand side of this product.
     */
    const T1& lhs() const { return this->m_lhs; }

    /**
     *  @return The right-hand side of this product.
     */
    const T2& rhs() const { return this->m_rhs; }

    /**
     *  Evaluate this function product.
     * 
     *  @param in           The argument for which the function is to be evaluated.
     *
     *  @return The function value of this function product for the given argument.
     */
    OutputType operator()(const InputType& in) const override {
        return this->m_lhs(in) * this->m_rhs(in);
    }
};


/*
 *  MARK: `Function` implementation
 */

/**
 *  A type that represents a mathematical function through its `operator()`.
 *
 *  @tparam _OutputType         The return type of the `operator()`.
 *  @tparam _InputType          The input type of the `operator()`.
 */
template <typename _OutputType, typename _InputType>
class Function {
public:
    // The return type of the `operator()`.
    using OutputType = _OutputType;

    // The input type of the `operator()`.
    using InputType = _InputType;


public:
    /*
     *  MARK: Destructor
     */

    /**
     *  The default destructor.
     */
    virtual ~Function() = default;


    /*
     *  MARK: `Function` behavior
     */

    /**
     *  Evaluate this function for the given argument.
     * 
     *  @param in        The argument at which the function is to be evaluated.
     *
     *  @return The function value for the given argument.
     */
    virtual OutputType operator()(const InputType& in) const = 0;
};


/*
 *  MARK: Aliases related to `Function`.
 */

/**
 *  A SFINAE expression that checks if the type T derives from `Function`.
 */
template <typename T>
using IsFunction = enable_if_t<std::is_base_of<Function<typename T::OutputType, typename T::InputType>, T>::value>;


/**
 *  Multiply one function by another.
 * 
 *  @tparam F1              The type of the first function.
 *  @tparam F2              The type of the second function.
 * 
 *  @note This function is only enabled for actual functions, i.e. functions that derive from `Function`.
 */
template <typename F1, typename F2, typename = IsFunction<F1>, typename = IsFunction<F2>>
FunctionProduct<F1, F2> operator*(const F1& lhs, const F2& rhs) {

    return FunctionProduct<F1, F2>(lhs, rhs);
}


}  // namespace GQCP
