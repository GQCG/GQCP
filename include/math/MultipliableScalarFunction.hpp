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
#ifndef MultipliableScalarFunction_hpp
#define MultipliableScalarFunction_hpp


#include "math/ScalarFunctionProduct.hpp"


namespace GQCP {


/**
 *  An extension of a scalar function that implements its operator* most generally as returning a ScalarFunctionProduct, whose operator() is implemented as the product of the underlying operator()s
 */
template <typename _Valued, typename _Scalar, int _Cols>
class MultipliableScalarFunction : public ScalarFunction<_Valued, _Scalar, _Cols> {
public:

    // Make sure the template parameters can be accessed like MultipliableScalarFunction::Valued
    using Valued = _Valued;
    using Scalar = _Scalar;
    enum {
        Cols = _Cols,
    };


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


#endif  /* MultipliableScalarFunction_h */
