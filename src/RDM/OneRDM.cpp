#include "RDM/OneRDM.hpp"


namespace GQCG {



/*
 *  CONSTRUCTORS
 */

OneRDM::OneRDM(Eigen::MatrixXd one_rdm) :
    BaseRDM (one_rdm.cols()),
    one_rdm (one_rdm)
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

/**
 *  diagonalises this->one_rdm and @returns the eigenvectors
 */
Eigen::MatrixXd OneRDM::diagonalise() {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes (this->one_rdm);
    this->one_rdm = saes.eigenvalues().asDiagonal();
    return saes.eigenvectors();
}


}  // namespace GQCG
