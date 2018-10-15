#include "RDM/TwoRDM.hpp"


namespace GQCG {


/*
 *  CONSTRUCTORS
 */

TwoRDM::TwoRDM(Eigen::Tensor<double, 4> two_rdm) :
    BaseRDM (two_rdm.dimensions()[0]),
    two_rdm (two_rdm)
{}

/**
 * Constructor where two_rdm = @param two_rdm_aaaa + @param two_rdm_bbbb + @param two_rdm_aabb  + @param two_rdm_bbaa
 */
TwoRDM::TwoRDM(Eigen::Tensor<double, 4> two_rdm_aaaa, Eigen::Tensor<double, 4> two_rdm_bbbb, Eigen::Tensor<double, 4> two_rdm_aabb, Eigen::Tensor<double, 4> two_rdm_bbaa) :
    BaseRDM (two_rdm_aaaa.dimensions()[0]),
    two_rdm (two_rdm_aaaa + two_rdm_bbbb + two_rdm_aabb + two_rdm_bbaa),
    two_rdm_aaaa (two_rdm_aaaa),
    two_rdm_bbbb (two_rdm_bbbb),
    two_rdm_aabb (two_rdm_aabb),
    two_rdm_bbaa (two_rdm_bbaa)
{}

}  // namespace GQCG

