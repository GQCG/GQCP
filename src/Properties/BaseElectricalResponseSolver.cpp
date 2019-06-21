#include "Properties/BaseElectricalResponseSolver.hpp"

#include "Eigen/Householder"


namespace GQCP {



/*
 *  PUBLIC METHODS
 */

/**
 *  Solve the linear response equations for the wave function response
 * 
 *  @param ham_par                  the Hamiltonian parameters
 *  @param dipole_integrals         the dipole integrals in an orthonormal orbital basis
 * 
 *  @return the wave function response
 */
VectorX<double> BaseElectricalResponseSolver::calculateWaveFunctionResponse(const HamiltonianParameters<double>& ham_par, const std::array<OneElectronOperator<double>, 3>& dipole_integrals) {

    const auto k_p = this->calculateParameterResponseConstant(ham_par);  // p for parameter
    const auto F_p = this->calculateParameterResponseForce(dipole_integrals);  // has 3 columns


    // This function is basically a wrapper around solving k_p x = -F_p
    Eigen::HouseholderQR<Eigen::MatrixXd> linear_solver (k_p);
    const VectorX<double> x = linear_solver.solve(-F_p);

    return x;
}


}  // namespace GQCP
