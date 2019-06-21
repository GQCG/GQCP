#include "Properties/CIElectricalResponseSolver.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param wf           the CI wave function
 *  @param molecule                 the molecule
 */
CIElectricalResponseSolver::CIElectricalResponseSolver(const WaveFunction& wave_function, const Molecule& molecule) :
    wave_function (wave_function),
    molecule (molecule)
{}



/*
 *  PUBLIC OVERRIDDEN METHODS
 */

/**
 *  @param ham_par                  the Hamiltonian parameters
 * 
 *  @return the parameter response constant (k_p), i.e. the second-order parameter partial derivative of the CI energy function
 */
SquareMatrix<double> CIElectricalResponseSolver::calculateParameterResponseConstant(const HamiltonianParameters<double>& ham_par) const {

    // k_p for DOCI is just the Hamiltonian evaluated in the Fock space
    return 2 * this->wave_function.get_fock_space().evaluateOperatorDense(ham_par, true);  // true: need to calculate diagonal values as well
}


/**
 *  @param dipole_integrals         the dipole integrals in an orthonormal orbital basis
 * 
 *  @return the parameter response force (F_p), i.e. the first-order parameter partial derivative of the perturbation derivative of the CI energy function
 */
Matrix<double, Dynamic, 3> CIElectricalResponseSolver::calculateParameterResponseForce(const std::array<OneElectronOperator<double>, 3>& dipole_integrals) const {

    const auto dim = this->wave_function.get_fock_space().get_dimension();  // dimension of the Fock space
    Matrix<double, Dynamic, 3> F_p = Matrix<double, Dynamic, 3>::Zero(dim, 3);


    // F_p for CI is a matrix-vector product of the total dipole operator and the CI wave function
    const auto nuclear_dipole_moment = this->molecule.calculateNuclearDipoleMoment();
    for (size_t m = 0; m < 3; m++) {  // m loops over the components of the electrical dipole

        // Generate the total dipole operator matrix representation
        auto mu_m = this->wave_function.get_fock_space().evaluateOperatorDense(dipole_integrals[m], true);  // true: need to calculate diagonal values as well

        // Set the nuclear contributions on the diagonal
        for (size_t J = 0; J < dim; J++) {
            mu_m(J,J) += nuclear_dipole_moment(m);
        }

        // Do the matvec explicitly and set the result as a column
        F_p.col(m) = mu_m * this->wave_function.get_coefficients();
    }


    return -2 * F_p;
}


}  // namespace GQCP
