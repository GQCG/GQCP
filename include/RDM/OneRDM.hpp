#ifndef GQCP_ONERDM_HPP
#define GQCP_ONERDM_HPP


#include "RDM/BaseRDM.hpp"

#include <Eigen/Dense>


namespace GQCP {

/**
 *  A class that holds the matrix representation of a 1-RDM
 */
class OneRDM : public BaseRDM {
private:
    Eigen::MatrixXd D;


public:
    // CONSTRUCTORS
    explicit OneRDM(const Eigen::MatrixXd& D);


    // GETTERS
    Eigen::MatrixXd get_matrix_representation() const { return this->D; }
    double get(size_t p, size_t q) const { return this->D(p, q); }


    // PUBLIC METHODS
    /**
     *  @return the 1-RDM's trace
     */
    double trace();

    /**
     *  diagonalizes the 1-RDM and @returns the eigenvectors
     */
    Eigen::MatrixXd diagonalize();
};


}  // namespace GQCP


#endif  // GQCP_ONERDM_HPP
