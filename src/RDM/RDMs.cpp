#include "RDM/RDMs.hpp"


namespace GQCG {


/*
 *  CONSTRUCTORS
 */

/**
 *  Constructor with @param one_rdm
 *  where one_rdm_aa and one_rdm_bb are @param one_rdm/2
 */
OneRDMs::OneRDMs(const OneRDM& one_rdm) :
        one_rdm (one_rdm),
        one_rdm_aa (one_rdm.get_matrix_representation()/2),
        one_rdm_bb (one_rdm.get_matrix_representation()/2)
{}

/**
 *  Constructor with @param one_rdm_aa and @param one_rdm_bb
 *  were one_rdm = @param one_rdm_aa + @param one_rdm_bb
 */
OneRDMs::OneRDMs(const OneRDM& one_rdm_aa, const OneRDM& one_rdm_bb) :
        one_rdm (one_rdm_aa.get_matrix_representation() + one_rdm_bb.get_matrix_representation()),
        one_rdm_aa (one_rdm_aa),
        one_rdm_bb (one_rdm_bb)
{}

/**
 *  Constructor
 *  where two_rdm = @param two_rdm_aaaa + @param two_rdm_aabb + @param two_rdm_bbaa  + @param two_rdm_bbbb
 */
TwoRDMs::TwoRDMs(const TwoRDM& two_rdm_aaaa, const TwoRDM& two_rdm_aabb, const TwoRDM& two_rdm_bbaa, const TwoRDM& two_rdm_bbbb) :
        two_rdm (two_rdm_aaaa.get_matrix_representation() + two_rdm_aabb.get_matrix_representation() + two_rdm_bbaa.get_matrix_representation() + two_rdm_bbbb.get_matrix_representation()),
        two_rdm_aaaa (two_rdm_aaaa),
        two_rdm_aabb (two_rdm_aabb),
        two_rdm_bbaa (two_rdm_bbaa),
        two_rdm_bbbb (two_rdm_bbbb)
{}


}  // namespace GQCG