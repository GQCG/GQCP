#include "geminals/AP1roGBivariationalSolver.hpp"

#include "geminals/AP1roGPSESolver.hpp"


namespace GQCP {


/*
 * CONSTRUCTORS
 */

/**
 *  @param N_P          the number of electrons
 *  @param ham_par      Hamiltonian parameters in an orthonormal orbital basis
 *  @param G            the initial guess for the AP1roG gemial coefficients
 *  @param extra_eq     the specification of the extra equation
 */
AP1roGBivariationalSolver::AP1roGBivariationalSolver(size_t N_P, const HamiltonianParameters& ham_par, const AP1roGGeminalCoefficients& G, ExtraEquation extra_eq) :
    BaseAP1roGSolver(N_P, ham_par, G),
    extra_eq (extra_eq)
{}


/**
 *  @param N_P          the number of electrons
 *  @param ham_par      Hamiltonian parameters in an orthonormal orbital basis
 *  @param extra_eq     the specification of the extra equation
 *
 *  The initial guess for the geminal coefficients is zero
 */
AP1roGBivariationalSolver::AP1roGBivariationalSolver(size_t N_P, const HamiltonianParameters& ham_par, ExtraEquation extra_eq) :
    BaseAP1roGSolver(N_P, ham_par),
    extra_eq (extra_eq)
{}


/**
 *  @param molecule     the molecule used for the AP1roG calculation
 *  @param ham_par      Hamiltonian parameters in an orthonormal orbital basis
 *  @param G            the initial guess for the AP1roG gemial coefficients
 *  @param extra_eq     the specification of the extra equation
 */
AP1roGBivariationalSolver::AP1roGBivariationalSolver(const Molecule& molecule, const HamiltonianParameters& ham_par, const AP1roGGeminalCoefficients& G, ExtraEquation extra_eq) :
    BaseAP1roGSolver(molecule, ham_par, G),
    extra_eq (extra_eq)
{}


/**
 *  @param molecule     the molecule used for the AP1roG calculation
 *  @param ham_par      Hamiltonian parameters in an orthonormal orbital basis
 *  @param extra_eq     the specification of the extra equation
 *
 *  The initial guess for the geminal coefficients is zero
 */
AP1roGBivariationalSolver::AP1roGBivariationalSolver(const Molecule& molecule, const HamiltonianParameters& ham_par, ExtraEquation extra_eq) :
    BaseAP1roGSolver(molecule, ham_par),
    extra_eq (extra_eq)
{}



/*
 *  PUBLIC METHODS
 */
void AP1roGBivariationalSolver::solve() {

    auto g = this->ham_par.get_g();


    // Solve the PSEs and set part of the solutions
    AP1roGPSESolver pse_solver (this->N_P, this->ham_par, this->geminal_coefficients);
    pse_solver.solve();

    this->geminal_coefficients = pse_solver.get_geminal_coefficients();
    this->electronic_energy = pse_solver.get_electronic_energy();


    // Initialize and solve the linear system Aq=b
    size_t dim = 1 + this->N_P * (this->K - N_P);

    Eigen::VectorXd b = Eigen::VectorXd::Zero(dim);
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(dim, dim);

    //      Initialize the extra equation
    switch (this->extra_eq) {
        case ExtraEquation::q0: {
            A(0,0) = 1;
            b(0) = 1;
            break;
        }

        case ExtraEquation::norm: {
            A(0,0) = 1;

            auto G = this->geminal_coefficients;
            for (size_t j = 0; j < this->N_P; j++) {
                for (size_t b = this->N_P; b < this->K; b++) {
                    size_t column_vector_index = G.vectorIndex(j, b);

                    A(0, 1 + column_vector_index) = G(j, b);
                }
            }

            b(0) = 1;
            break;
        }

        default:
            break;
    }

    //      Initialize the rest of A
    auto J = pse_solver.calculateJacobian(this->geminal_coefficients.asVector());

    for (size_t i = 0; i < this->N_P; i++) {
        for (size_t a = this->N_P; a < this->K; a++) {
            size_t row_vector_index = this->geminal_coefficients.vectorIndex(i, a);

            // First column
            A(1 + row_vector_index, 0) = g(i,a,i,a);

            // Large lower right block
            for (size_t j = 0; j < this->N_P; j++) {
                for (size_t b = this->N_P; b < this->K; b++) {
                    size_t column_vector_index = this->geminal_coefficients.vectorIndex(j, b);

                    A(1 + row_vector_index, 1 + column_vector_index) = J(column_vector_index, row_vector_index) + g(i,a,i,a) * this->geminal_coefficients(j, b);  // transpose of the Jacobian
                }
            }  // j and b

        }
    }

    Eigen::HouseholderQR<Eigen::MatrixXd> linear_solver (A);
    Eigen::VectorXd q = linear_solver.solve(b);
    assert(std::abs((A * q).norm() - 1) < 1.0e-12);


    // Set the solution
    this->bivariational_coefficients.q0 = q(0);
    this->bivariational_coefficients.q = AP1roGVariables(q.tail(this->N_P * (this->K - N_P)), this->N_P, this->K);
}


}  // namespace GQCP
