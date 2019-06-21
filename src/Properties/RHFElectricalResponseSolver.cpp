#include "Properties/RHFElectricalResponseSolver.hpp"

#include "RHF/RHF.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param N_P          the number of electron pairs
 */
RHFElectricalResponseSolver::RHFElectricalResponseSolver(const size_t N_P) :
    N_P (N_P)
{}



/*
 *  PUBLIC OVERRIDDEN METHODS
 */

/**
 *  @param ham_par                  the Hamiltonian parameters
 * 
 *  @return the parameter response constant (k_p), i.e. the second-order parameter partial derivative of the RHF energy function
 */
SquareMatrix<double> RHFElectricalResponseSolver::calculateParameterResponseConstant(const HamiltonianParameters<double>& ham_par) const {

    // k_p for RHF is the RHF orbital Hessian
    const auto RHF_orbital_hessian_tensor = calculateRHFOrbitalHessianTensor(ham_par, this->N_P);

    // We have now calculated the full KxKxKxK tensor (K: number of spatial orbitals)
    // We should reduce the tensor dimension to only include the non-redundant rotations (ignore zero elements) and transform to a Hessian matrix

    return RHF_orbital_hessian_tensor.pairWiseReduce(this->N_P, 0, this->N_P, 0);  // compound both indices (ai, bj) into (m, n), knowing that virtual indices start at N_P
}


/**
 *  @param dipole_integrals         the dipole integrals in an orthonormal orbital basis
 * 
 *  @return the parameter response force (F_p), i.e. the first-order parameter partial derivative of the perturbation derivative of the RHF energy function
 */
Matrix<double, Dynamic, 3> RHFElectricalResponseSolver::calculateParameterResponseForce(const std::array<OneElectronOperator<double>, 3>& dipole_integrals) const {

    const auto K = dipole_integrals[0].get_K();


    const auto dim = K * this->N_P;  // number of non-redundant orbital rotation generators
    Matrix<double, Dynamic, 3> F_p = Matrix<double, Dynamic, 3>::Zero(dim, 3);


    for (size_t m = 0; m < 3; m++) {  // m loops over the components of the electrical dipole

        // mu_m is a KxK matrix representation of the m-th component of the electrical electron dipole moment
        // Put the virtual-occupied mu_m(a,i) elements into the response force vector
        F_p.col(m) = dipole_integrals[m].pairWiseReduce(this->N_P, 0);
    }

    return 4 * F_p;
}


}  // namespace GQCP
