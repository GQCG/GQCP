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
#ifndef ScalarFunctionProduct_hpp
#define ScalarFunctionProduct_hpp


#include <type_traits>

#include "math/ScalarFunction.hpp"


namespace GQCP {


/**
 *  A class template that represents a product of scalar functions (of the same type)
 *
 *  operator() returns the product of the operator()s of both or the arguments
 *
 *  @tparam T    the types of the scalar functions
 */
template <typename T>
class ScalarFunctionProduct : public ScalarFunction<typename T::Valued, typename T::Scalar, T::Cols> {
    static_assert(std::is_base_of<ScalarFunction<typename T::Valued, typename T::Scalar, T::Cols>, T>::value, "ScalarFunctionProduct: T must inherit from ScalarFunction");


private:
    const T& lhs;
    const T& rhs;


public:

    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param lhs      the left-hand side of the product
     *  @param rhs      the right-hand side of the product
     */
    ScalarFunctionProduct(const T& lhs, const T& rhs) :
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
    typename T::Valued operator()(const Eigen::Matrix<typename T::Scalar, T::Cols, 1>& x) const override {
        return this->lhs(x) * this->rhs(x);
    }
};


}  // namespace GQCP


#endif  /* ScalarFunctionProduct_h */
