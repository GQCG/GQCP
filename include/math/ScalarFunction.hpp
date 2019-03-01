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


#include "typedefs.hpp"


namespace GQCP {


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

    // Make sure the template parameters can be accessed like ScalarFunction::Valued
    using Valued = _Valued;
    using Scalar = _Scalar;
    enum {
        Cols = _Cols,
    };


    /**
     *  @param x        the vector/point at which the scalar function is to be evaluated
     *
     *  @return the scalar function value at the given point
     */
    virtual _Valued operator()(const Vector<_Scalar, _Cols>& x) const = 0;
};


}  // namespace GQCP



#endif  /* ScalarFunction_hpp */
