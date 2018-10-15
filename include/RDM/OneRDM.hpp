#ifndef GQCG_ONERDM_HPP
#define GQCG_ONERDM_HPP


#include "RDM/BaseRDM.hpp"

#include <Eigen/Dense>


namespace GQCG {


class OneRDM : public BaseRDM {
private:
    Eigen::MatrixXd oneRDM;


public:
    // CONSTRUCTORS
    OneRDM(Eigen::MatrixXd oneRDM);


    // GETTERS
    Eigen::MatrixXd get_matrix_representation() const { return this->oneRDM; }
    double get(size_t p, size_t q) const { return this->oneRDM(p, q); }
};


}  // namespace GQCG


#endif  // GQCG_ONERDM_HPP
