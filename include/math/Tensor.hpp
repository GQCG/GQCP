#ifndef Tensor_h
#define Tensor_h


#include <unsupported/Eigen/CXX11/Tensor>


namespace GQCP {


/**
 *  An extension of the Eigen::Tensor class, with extra operations
 *
 *  We have decided to inherit from Eigen::Tensor, because we will use different hierarchies: see also: https://eigen.tuxfamily.org/dox-devel/TopicCustomizing_InheritingMatrix.html
 */
template <typename _Scalar, int _Rank>
class Tensor : public Eigen::Tensor<_Scalar, _Rank> {
public:

    using Scalar = _Scalar;
    enum {
        Rank = _Rank
    };

    using Self = Tensor<Scalar, Rank>;
    using Base = Eigen::Tensor<Scalar, Rank>;


public:

    /*
     *  CONSTRUCTORS
     */
    using Base::Base;  // inherit Base constructors
};


}  // namespace GQCP


#endif  /* Tensor_h */
