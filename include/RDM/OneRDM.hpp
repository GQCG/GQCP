#ifndef GQCG_ONERDM_HPP
#define GQCG_ONERDM_HPP


#include "RDM/BaseRDM.hpp"

#include <Eigen/Dense>
#include <Eigenpair.hpp>


namespace GQCG {

/**
 *  A class that holds the matrix representations of a 1RDM
 */
class OneRDM : public BaseRDM {
private:
    Eigen::MatrixXd one_rdm;


public:
    // CONSTRUCTORS
    explicit OneRDM(Eigen::MatrixXd one_rdm);


    // GETTERS
    Eigen::MatrixXd get_matrix_representation() const { return this->one_rdm; }
    double get(size_t p, size_t q) const { return this->one_rdm(p, q); }


    // PUBLIC METHODS
    /**
     *  @return the trace of this->one_rdm
     */
    double trace();

    /**
     *  diagonalises this->one_rdm and @returns the eigenvectors
     */
    Eigen::MatrixXd diagonalise();
};


}  // namespace GQCG


#endif  // GQCG_ONERDM_HPP
