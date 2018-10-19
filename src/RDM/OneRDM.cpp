#include "RDM/OneRDM.hpp"


namespace GQCP {



/*
 *  CONSTRUCTORS
 */

OneRDM::OneRDM(const Eigen::MatrixXd& D) :
    BaseRDM (D.cols()),
    D (D)
{
    // Check if the 1-RDM is represented as a square matrix
    if (D.cols() != D.rows()) {
        throw std::invalid_argument("1-RDMs have to be represented as a square matrix.");
    }
}



/*
 *  PUBLIC METHODS
 */

/**
 *  @return the 1-RDM's trace
 */
double OneRDM::trace(){
    return this->D.trace();
}


/**
 *  diagonalizes the 1-RDM and @returns the eigenvectors
 */
Eigen::MatrixXd OneRDM::diagonalize() {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes (this->D);
    this->D = saes.eigenvalues().asDiagonal();
    return saes.eigenvectors();
}


}  // namespace GQCP
