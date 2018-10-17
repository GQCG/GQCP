#ifndef GQCG_TWORDM_HPP
#define GQCG_TWORDM_HPP


#include "RDM/BaseRDM.hpp"

#include <unsupported/Eigen/CXX11/Tensor>


namespace GQCG {

/**
 *  A class that holds the tensor representations of a 2RDM
 */
class TwoRDM : public BaseRDM {
private:
    Eigen::Tensor<double, 4> d;


public:
    // CONSTRUCTORS
    explicit TwoRDM(const Eigen::Tensor<double, 4>& d);


    // GETTERS
    Eigen::Tensor<double, 4> get_matrix_representation() const { return this->d; }
    double get(size_t p, size_t q, size_t r, size_t s) const { return this->d(p, q, r, s); }


    // PUBLIC METHODS
    /**
     *  @return the trace of the 2-RDM @param: d(p,p,q,q)
     */
    double trace();

    /**
     *  @return a partial contraction of the 2-RDM,
     *  where D(p,q) = d(p,q,r,r)
     */
    Eigen::MatrixXd reduce();
};


}  // namespace GQCG


#endif  // GQCG_TWORDM_HPP
