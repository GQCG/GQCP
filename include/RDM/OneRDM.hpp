#ifndef GQCG_ONERDM_HPP
#define GQCG_ONERDM_HPP


#include "RDM/BaseRDM.hpp"

#include <Eigen/Dense>


namespace GQCG {

/**
 *  A class that holds the matrix representations of a 1RDM
 */
class OneRDM : public BaseRDM {
private:
    Eigen::MatrixXd one_rdm;  // spin-summed (total) 1-RDM

    Eigen::MatrixXd one_rdm_aa;  // alpha-alpha (a-a) 1-RDM
    Eigen::MatrixXd one_rdm_bb;  // beta-beta (b-b) 1-RDM


public:
    // CONSTRUCTORS
    /**
     * Constructor with @param one_rdm where one_rdm_aa and one_rdm_bb are @param one_rdm/2
     */
    OneRDM(Eigen::MatrixXd one_rdm);

    /**
     * Constructor with @param one_rdm_aa and @param one_rdm_bb were one_rdm = @param one_rdm_aa + @param one_rdm_bb
     */
    OneRDM(Eigen::MatrixXd one_rdm_aa, Eigen::MatrixXd one_rdm_bb);


    // GETTERS
    Eigen::MatrixXd get_matrix_representation() const { return this->one_rdm; }
    double get(size_t p, size_t q) const { return this->one_rdm(p, q); }


    // PUBLIC METHODS
    /**
     *  @return the trace of this->one_rdm
     */
    double trace();
};


}  // namespace GQCG


#endif  // GQCG_ONERDM_HPP
