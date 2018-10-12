#ifndef GQCG_WAVEFUNCTION_HPP
#define GQCG_WAVEFUNCTION_HPP


#include "FockSpace/BaseFockSpace.hpp"

#include "common.hpp"



namespace GQCG {

using BaseFockSpace_sptr = std::shared_ptr<BaseFockSpace>;

class WaveFunction {
private:
    BaseFockSpace_sptr base_fock_space;
    Eigen::MatrixXd eigenvectors;
    Eigen::VectorXd eigenvalues;

public:
    WaveFunction(BaseFockSpace& base_fock_space, Eigen::VectorXd eigenvector, double eigenvalue);
    WaveFunction(BaseFockSpace& base_fock_space, Eigen::MatrixXd eigenvector, Eigen::VectorXd eigenvalue);

    // GETTERS
    double get_eigenvalue(size_t index = 0) { return eigenvalues(index); }
    Eigen::VectorXd get_eigenvalues() { return eigenvalues; }
    Eigen::VectorXd get_eigenvector(size_t index = 0) { return eigenvectors.col(index); }
    Eigen::MatrixXd get_eigenvectors() { return eigenvectors; }
    BaseFockSpace& get_fock_space() { return *base_fock_space; }
};


}  // namespace GQCG
#endif  // GQCG_WAVEFUNCTION_HPP
