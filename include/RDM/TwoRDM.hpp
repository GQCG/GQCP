#ifndef GQCG_TWORDM_HPP
#define GQCG_TWORDM_HPP


#include "RDM/BaseRDM.hpp"

#include <unsupported/Eigen/CXX11/Tensor>


namespace GQCG {


class TwoRDM : public BaseRDM {
private:
    Eigen::Tensor<double, 4> two_rdm;

    Eigen::Tensor<double, 4> two_rdm_aaaa;
    Eigen::Tensor<double, 4> two_rdm_aabb;
    Eigen::Tensor<double, 4> two_rdm_bbaa;
    Eigen::Tensor<double, 4> two_rdm_bbbb;


public:
    // CONSTRUCTORS
    TwoRDM(Eigen::Tensor<double, 4> two_rdm);

    /**
     * Constructor where two_rdm = @param two_rdm_aaaa + @param two_rdm_bbbb + @param two_rdm_aabb  + @param two_rdm_bbaa
     */
    TwoRDM(Eigen::Tensor<double, 4> two_rdm_aaaa, Eigen::Tensor<double, 4> two_rdm_bbbb, Eigen::Tensor<double, 4> two_rdm_aabb, Eigen::Tensor<double, 4> two_rdm_bbaa);


    // GETTERS
    Eigen::Tensor<double, 4> get_tensor_representation() const { return this->two_rdm; }
    double get(size_t p, size_t q, size_t r, size_t s) const { return this->two_rdm(p, q, r, s); }
};


}  // namespace GQCG


#endif  // GQCG_TWORDM_HPP
