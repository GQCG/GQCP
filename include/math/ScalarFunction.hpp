#ifndef ScalarFunction_hpp
#define ScalarFunction_hpp


#include "typedefs.hpp"


namespace GQCP {


/**
 *  A class template representing a mathematical scalar function through overriding operator()
 *
 *  @tparam _Valued     the return type of the scalar function
 *  @tparam _Scalar     the type of the scalars of the input vector
 *  @tparam _Cols       the dimension of the input vector, an integer, or Eigen::Dynamic representing an unknown number of columns at compile time
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
     *  operator() implements the notion of a `_Valued`-valued scalar function, accepting a vector of `_Scalar`s
     */
    virtual _Valued operator()(const Vector<_Scalar, _Cols>& x) const = 0;
};


}  // namespace GQCP



#endif /* ScalarFunction_hpp */
