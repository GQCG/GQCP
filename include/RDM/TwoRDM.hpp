#ifndef GQCG_TWORDM_HPP
#define GQCG_TWORDM_HPP


#include "RDM/BaseRDM.hpp"

#include <unsupported/Eigen/CXX11/Tensor>


namespace GQCG {


class TwoRDM : public BaseRDM {
private:
    Eigen::Tensor<double, 4> twoRDM;


public:
    // CONSTRUCTORS
    TwoRDM(Eigen::Tensor<double, 4> twoRDM);


    // GETTERS
    Eigen::Tensor<double, 4> get_tensor_representation() const { return this->twoRDM; }
    double get(size_t p, size_t q, size_t r, size_t s) const { return this->twoRDM(p, q, r, s); }
};


}  // namespace GQCG


#endif  // GQCG_TWORDM_HPP
