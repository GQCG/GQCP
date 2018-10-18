#include "DOCINewtonOrbitalOptimizer.hpp"


namespace GQCG {


/*
 *  CONSTRUCTORS
 */
/**
 *  Constructor based on a given @param fock_space, Hamiltonian parameters @param ham_par, @param solver_options for solving the DOCI eigenvalue problem, a @param oo_convergence_threshold and a @param maximum_number_of_oo_iterations
 */
DOCINewtonOrbitalOptimizer::DOCINewtonOrbitalOptimizer(const GQCG::FockSpace& fock_space, const GQCG::HamiltonianParameters& ham_par, numopt::eigenproblem::BaseSolverOptions& solver_options, double oo_convergence_threshold, size_t maximum_number_of_oo_iterations) :
    fock_space (fock_space),
    ham_par (ham_par),
    solver_options (solver_options),
    oo_convergence_threshold (oo_convergence_threshold),
    maximum_number_of_oo_iterations (maximum_number_of_oo_iterations)
{}



/*
 *  PUBLIC METHODS
 */
void DOCINewtonOrbitalOptimizer::solve() {

    bool is_OO_converged = false;
    size_t OO_iterations = 0;
    while (!is_OO_converged) {

        // Solve the DOCI eigenvalue equation, using the options provided
        this->solve(solver_options_ptr);


        // Calculate the 1- and 2-RDMs
        this->calculate1RDMs();
        this->calculate2RDMs();


        // Calculate the electronic gradient at kappa = 0
        Eigen::MatrixXd F = this->so_basis.calculateGeneralizedFockMatrix(this->one_rdm, this->two_rdm);
        Eigen::MatrixXd gradient_matrix = 2 * (F - F.transpose());
        Eigen::VectorXd gradient_vector = cpputil::linalg::strictLowerTriangle(gradient_matrix);  // gradient vector with the free parameters, at kappa = 0


        // Calculate the electronic Hessian at kappa = 0
        Eigen::Tensor<double, 4> W = this->so_basis.calculateSuperGeneralizedFockMatrix(this->one_rdm, this->two_rdm);
        Eigen::Tensor<double, 4> hessian_tensor (this->K, this->K, this->K, this->K);
        hessian_tensor.setZero();

        for (size_t p = 0; p < this->K; p++) {
            for (size_t q = 0; q < this->K; q++) {
                for (size_t r = 0; r < this->K; r++) {
                    for (size_t s = 0; s < this->K; s++) {
                        hessian_tensor(p,q,r,s) = W(p,q,r,s) - W(p,q,s,r) + W(q,p,s,r) - W(q,p,r,s) + W(r,s,p,q) - W(r,s,q,p) + W(s,r,q,p) - W(s,r,p,q);
                    }
                }
            }
        }
        Eigen::MatrixXd hessian_matrix = cpputil::linalg::strictLowerTriangle(hessian_tensor);  // hessian matrix with only the free parameters, at kappa = 0

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> hessian_solver (hessian_matrix);


        // At this moment, we have calculated the electronic gradient and electronic Hessian at kappa = 0
        // Perform a Newton-step to find orbital rotation parameters kappa
        numopt::GradientFunction gradient_function = [gradient_vector](const Eigen::VectorXd& x) { return gradient_vector; };
        numopt::JacobianFunction hessian_function = [hessian_matrix](const Eigen::VectorXd& x) { return hessian_matrix; };

        Eigen::VectorXd kappa_vector = numopt::newtonStep(Eigen::VectorXd::Zero(this->K), gradient_function, hessian_function);  // with only the free parameters


        // If the calculated norm is zero, we have reached a critical point
        if (gradient_vector.norm() < OO_convergence_threshold) {

            // If we have found a critical point, but we have a negative eigenvalue for the Hessian, continue in that direction
            if (hessian_solver.eigenvalues()(0) < 0) {
                kappa_vector = hessian_solver.eigenvectors().col(0);
            }
            else {  // the Hessian is confirmed to be positive definite, so we have reached a minimum
                is_OO_converged = true;
            }


        } else {
            OO_iterations++;

            if (OO_iterations >= maximum_number_of_OO_iterations) {
                throw std::runtime_error("DOCI::orbitalOptimize(): The OO-DOCI procedure failed to converge in the maximum number of allowed iterations.");
            }
        }


        // Change kappa back to a matrix
        Eigen::MatrixXd kappa_matrix = cpputil::linalg::fillStrictLowerTriangle(kappa_vector);  // containing all parameters, so this is in anti-Hermitian (anti-symmetric) form
        Eigen::MatrixXd kappa_matrix_transpose = kappa_matrix.transpose();  // store the transpose in an auxiliary variable to avoid aliasing issues
        kappa_matrix -= kappa_matrix_transpose;  // fillStrictLowerTriangle only returns the lower triangle, so we must construct the anti-Hermitian (anti-symmetric) matrix


        // Calculate the unitary rotation matrix that we can use to rotate the basis
        Eigen::MatrixXd U = (-kappa_matrix).exp();


        // Transform the integrals to the new orthonormal basis
        this->so_basis.rotate(U);  // this checks if U is actually unitary


        // If we're using a DavidsonSolver, we should update the the initial guesses to be the current eigenvectors
        if (solver_options_ptr->get_solver_type() == numopt::eigenproblem::SolverType::DAVIDSON) {
            auto davidson_solver_options_ptr = dynamic_cast<numopt::eigenproblem::DavidsonSolverOptions*>(solver_options_ptr);  // this now points to the used solver options

            for (size_t i = 0; i < davidson_solver_options_ptr->number_of_required_eigenpairs; i++) {
                davidson_solver_options_ptr->X_0.col(i) = this->get_eigenpair(i).get_eigenvector();
            }
        }
    }  // while not converged
}


}  // namespace GQCG
