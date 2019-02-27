#ifndef ScalarFunctionProduct_hpp
#define ScalarFunctionProduct_hpp


namespace GQCP {

/**
 *  @tparam T    the types of the scalar functions
 *
 *  TODO: make sure T1 and T2 are actually scalar functions (or derived classes)
 */
template <typename T>
class ScalarFunctionProduct : public ScalarFunction<typename T::Valued, typename T::Scalar, T::Cols> {
public:
    const T& lhs;
    const T& rhs;


    // CONSTRUCTORS
    ScalarFunctionProduct(const T& lhs, const T& rhs) :
        lhs (lhs),
        rhs (rhs)
    {}


    /**
     *  A product of scalar functions is again a scalar function, in which each argument should be evaluated at the given argument
     */
    typename T::Valued operator()(const Eigen::Matrix<typename T::Scalar, T::Cols, 1>& x) const override {
        return this->lhs(x) * this->rhs(x);
    }
};



}



#endif  /* ScalarFunctionProduct_h */
