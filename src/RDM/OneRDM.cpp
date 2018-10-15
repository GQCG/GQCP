#include "RDM/OneRDM.hpp"


namespace GQCG {



/*
 *  CONSTRUCTORS
 */

/**
 * Constructor with @param one_rdm where one_rdm_aa and one_rdm_bb are @param one_rdm/2
 */
OneRDM::OneRDM(Eigen::MatrixXd one_rdm) :
    BaseRDM (one_rdm.cols()),
    one_rdm (one_rdm),
    one_rdm_aa (one_rdm/2),
    one_rdm_bb (one_rdm/2)
{}

/**
 * Constructor with @param one_rdm_aa and @param one_rdm_bb were one_rdm = @param one_rdm_aa + @param one_rdm_bb
 */
OneRDM::OneRDM(Eigen::MatrixXd one_rdm_aa, Eigen::MatrixXd one_rdm_bb) :
    BaseRDM (one_rdm_aa.cols()),
    one_rdm (one_rdm_aa + one_rdm_bb),
    one_rdm_aa (one_rdm_aa),
    one_rdm_bb (one_rdm_bb)
{}

/*
 *  PUBLIC METHODS
 */

/**
 *  @return the trace of this->one_rdm
 */
double OneRDM::trace(){
    return this->one_rdm.trace();
}

}  // namespace GQCG
