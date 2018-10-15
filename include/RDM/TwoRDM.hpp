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
    Eigen::Tensor<double, 4> two_rdm;  // spin-summed (total) 2-RDM

    Eigen::Tensor<double, 4> two_rdm_aaaa;  // a-a-a-a 2-RDM
    Eigen::Tensor<double, 4> two_rdm_aabb;  // a-a-b-b 2-RDM
    Eigen::Tensor<double, 4> two_rdm_bbaa;  // b-a-a-b 2-RDM
    Eigen::Tensor<double, 4> two_rdm_bbbb;  // b-b-b-b 2-RDM


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


    // PUBLIC METHODS
    /**
     *  @return the trace of this->two_rdm
     */
    double trace();

    /**
     *  @return Eigen::MatrixXd D, the reduced-over 2-RDM : D(p,q) = this->two_rdm(p,q,r,r)
     */
    Eigen::MatrixXd reduce_2RDM();
};


}  // namespace GQCG


#endif  // GQCG_TWORDM_HPP
