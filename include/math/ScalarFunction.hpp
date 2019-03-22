// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2019  the GQCG developers
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
#ifndef ScalarFunction_hpp
#define ScalarFunction_hpp


/*
 *  This file contains the source code for both ScalarFunctionProduct and ScalarFunction
 */

#include "math/Matrix.hpp"


namespace GQCP {


/*
 *  Forward declaration of ScalarFunction to be used in ScalarFunctionProduct
 */
template <typename _Valued, typename _Scalar, int _Cols>
class ScalarFunction;



/*
 *  ScalarFunctionProduct
 */

template <typename T, typename U>
using product_t = decltype(std::declval<T>() * std::declval<U>());


/**
 *  A class template that represents a product of scalar functions (of the same type), such that the evaluation of a product is the product of the evaluations
 *
 *  @tparam T1      the left-hand side scalar function type
 *  @tparam T2      the right-hand side scalar function type
 */
template <typename T1, typename T2 = T1>
    class ScalarFunctionProduct : public ScalarFunction<product_t<typename T1::Valued, typename T2::Valued>, typename T1::Scalar, T1::Cols> {

    static_assert(std::is_base_of<ScalarFunction<typename T1::Valued, typename T1::Scalar, T1::Cols>, T1>::value, "ScalarFunctionProduct: T1 must inherit from ScalarFunction");
    static_assert(std::is_base_of<ScalarFunction<typename T2::Valued, typename T2::Scalar, T2::Cols>, T2>::value, "ScalarFunctionProduct: T2 must inherit from ScalarFunction");

    static_assert(std::is_same<typename T1::Scalar, typename T2::Scalar>::value, "ScalarFunctionProduct: T1 and T2 should have the same Scalar type");
    static_assert(T1::Cols == T2::Cols, "ScalarFunctionProduct: T1 and T2 should have the same Cols");


private:
    const T1& lhs;  // left-hand side
    const T2& rhs;  // right-hand side


public:
    using Valued = product_t<typename T1::Valued, typename T2::Valued>;
    using Scalar = typename T1::Scalar;  // equal to T2::Scalar
    enum {
        Cols = T1::Cols  // equal to T2::Cols
    };


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param lhs      the left-hand side of the product
     *  @param rhs      the right-hand side of the product
     */
    ScalarFunctionProduct(const T1& lhs, const T2& rhs) :
        lhs (lhs),
        rhs (rhs)
    {}


    /*
     *  OPERATORS
     */

    /**
     *  @param x        the vector/point at which the scalar function product is to be evaluated
     *
     *  @return the product of the evaluated left-hand and right-hand side scalar functions
     */
    Valued operator()(const Vector<Scalar, Cols>& x) const override {
        return this->lhs(x) * this->rhs(x);
    }
};




/*
 *  ScalarFunction
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
    enum {
        Cols = _Cols,
    };


public:

    /*
     *  OPERATORS
     */
    
    /**
     *  @param x        the vector/point at which the scalar function is to be evaluated
     *
     *  @return the scalar function value at the given point
     */
    virtual _Valued operator()(const Vector<_Scalar, _Cols>& x) const = 0;


    /**
     *  @param rhs      the right-hand side of the product
     *
     *  @return the product of this and the right-hand side scalar function
     */
    ScalarFunctionProduct<ScalarFunction<_Valued, _Scalar, _Cols>> operator*(const ScalarFunction<_Valued, _Scalar, _Cols>& rhs) const {
        return ScalarFunctionProduct<ScalarFunction<_Valued, _Scalar, _Cols>>(*this, rhs);
    }
};


}  // namespace GQCP



#endif  /* ScalarFunction_hpp */
