#ifndef GQCG_WAVEFUNCTION_HPP
#define GQCG_WAVEFUNCTION_HPP


#include "common.hpp"



namespace GQCG {

using BaseFockSpace_sptr = std::shared_ptr<BaseFockSpace>;

class WaveFunction {
private:
    BaseFockSpace_sptr base_fock_space;
    Eigen::MatrixXd eigenvector;
    Eigen::VectorXd eigenvalue;


public:
    WaveFunction(BaseFockSpace& base_fock_space, Eigen::VectorXd eigenvector, double eigenvalue)
    WaveFunction(BaseFockSpace& base_fock_space, Eigen::MatrixXd eigenvector, Eigen::VectorXd eigenvalue)
};


}  // namespace GQCG
#endif  // GQCG_WAVEFUNCTION_HPP
