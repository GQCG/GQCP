#include "DOCINewtonOrbitalOptimizer.hpp"

#include <cpputil.hpp>
#include <unsupported/Eigen/MatrixFunctions>


#include "CISolver/CISolver.hpp"
#include "RDM/DOCIRDMBuilder.hpp"


namespace GQCG {


/*
 *  CONSTRUCTORS
 */
/**
 *  Constructor based on a given @param fock_space, Hamiltonian parameters @param ham_par, @param solver_options for solving the DOCI eigenvalue problem, a @param oo_convergence_threshold and a @param maximum_number_of_oo_iterations
 */
DOCINewtonOrbitalOptimizer::DOCINewtonOrbitalOptimizer(const GQCG::FockSpace& fock_space, const GQCG::HamiltonianParameters& ham_par, numopt::eigenproblem::BaseSolverOptions& solver_options, double oo_convergence_threshold, size_t maximum_number_of_oo_iterations) :
    doci (GQCG::DOCI (fock_space)),
    fock_space (fock_space),
    ham_par (ham_par),
    solver_options (solver_options),
    doci_solver (GQCG::CISolver (this->doci, this->ham_par)),
    oo_convergence_threshold (oo_convergence_threshold),
    maximum_number_of_oo_iterations (maximum_number_of_oo_iterations)
{}


/*
 *  GETTERS
 */
std::vector<numopt::eigenproblem::Eigenpair> DOCINewtonOrbitalOptimizer::get_eigenpairs() const {
    if (this->is_converged) {
        return this->doci_solver.get_eigenpairs();
    } else {
        throw std::logic_error("You are trying to get eigenpairs but the orbital optimization hasn't converged (yet).");
    }
}

numopt::eigenproblem::Eigenpair DOCINewtonOrbitalOptimizer::get_eigenpair(size_t index) const {
    if (this->is_converged) {
        return this->doci_solver.get_eigenpair(index);
    } else {
        throw std::logic_error("You are trying to get eigenpairs but the orbital optimization hasn't converged (yet).");
    }
}



/*
 *  PUBLIC METHODS
 */
/**
 *  Perform the orbital optimization
 */
void DOCINewtonOrbitalOptimizer::solve() {

    auto K = this->fock_space.get_K();


    size_t oo_iterations = 0;
    while (!(this->is_converged)) {

        // Solve the DOCI eigenvalue equation, using the options provided
        this->doci_solver = GQCG::CISolver(this->doci, this->ham_par);  // update the CI solver with the rotated Hamiltonian parameters
        this->doci_solver.solve(this->solver_options);
        GQCG::WaveFunction ground_state = this->doci_solver.get_wavefunction();

        // Calculate the 1- and 2-RDMs
        GQCG::DOCIRDMBuilder rdm_builder (this->fock_space);
        auto one_rdm = rdm_builder.calculate1RDMs(ground_state.get_coefficients()).one_rdm;  // spin-summed 1-RDM
        auto two_rdm = rdm_builder.calculate2RDMs(ground_state.get_coefficients()).two_rdm;  // spin-summed 2-RDM


        // Calculate the electronic gradient at kappa = 0
        Eigen::MatrixXd F = this->ham_par.calculateGeneralizedFockMatrix(one_rdm, two_rdm).get_matrix_representation();
        Eigen::MatrixXd gradient_matrix = 2 * (F - F.transpose());
        Eigen::VectorXd gradient_vector = cpputil::linalg::strictLowerTriangle(gradient_matrix);  // gradient vector with the free parameters, at kappa = 0


        // Calculate the electronic Hessian at kappa = 0
        Eigen::Tensor<double, 4> W = this->ham_par.calculateSuperGeneralizedFockMatrix(one_rdm, two_rdm).get_matrix_representation();
        Eigen::Tensor<double, 4> hessian_tensor (K, K, K, K);
        hessian_tensor.setZero();

        for (size_t p = 0; p < K; p++) {
            for (size_t q = 0; q < K; q++) {
                for (size_t r = 0; r < K; r++) {
                    for (size_t s = 0; s < K; s++) {
                        hessian_tensor(p,q,r,s) = W(p,q,r,s) - W(p,q,s,r) + W(q,p,s,r) - W(q,p,r,s) + W(r,s,p,q) - W(r,s,q,p) + W(s,r,q,p) - W(s,r,p,q);
                    }
                }
            }
        }
        Eigen::MatrixXd hessian_matrix = cpputil::linalg::strictLowerTriangle(hessian_tensor);  // hessian matrix with only the free parameters, at kappa = 0

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> hessian_solver (hessian_matrix);


        // Perform a Newton-step to find orbital rotation parameters kappa
        numopt::GradientFunction gradient_function = [gradient_vector](const Eigen::VectorXd& x) { return gradient_vector; };
        numopt::JacobianFunction hessian_function = [hessian_matrix](const Eigen::VectorXd& x) { return hessian_matrix; };

        Eigen::VectorXd kappa_vector = numopt::newtonStep(Eigen::VectorXd::Zero(K), gradient_function, hessian_function);  // with only the free parameters


        // If the calculated norm is zero, we have reached a critical point
        if (gradient_vector.norm() < this->oo_convergence_threshold) {

            // If we have found a critical point, but we have a negative eigenvalue for the Hessian, continue in that direction
            if (hessian_solver.eigenvalues()(0) < 0) {
                kappa_vector = hessian_solver.eigenvectors().col(0);
            }
            else {  // the Hessian is confirmed to be positive definite, so we have reached a minimum
                this->is_converged = true;
            }


        } else {
            oo_iterations++;

            if (oo_iterations >= this->maximum_number_of_oo_iterations) {
                throw std::runtime_error("DOCINewtonOrbitalOptimizer.solve(): The OO-DOCI procedure failed to converge in the maximum number of allowed iterations.");
            }
        }


        // Change kappa back to a matrix
        Eigen::MatrixXd kappa_matrix = cpputil::linalg::fillStrictLowerTriangle(kappa_vector);  // containing all parameters, so this is in anti-Hermitian (anti-symmetric) form
        Eigen::MatrixXd kappa_matrix_transpose = kappa_matrix.transpose();  // store the transpose in an auxiliary variable to avoid aliasing issues
        kappa_matrix -= kappa_matrix_transpose;  // fillStrictLowerTriangle only returns the lower triangle, so we must construct the anti-Hermitian (anti-symmetric) matrix


        // Calculate the unitary rotation matrix that we can use to rotate the basis
        Eigen::MatrixXd U = (-kappa_matrix).exp();


        // Transform the integrals to the new orthonormal basis
        this->ham_par.rotate(U);  // this checks if U is actually unitary


        // If we're using a Davidson solver, we should update the initial guesses in the solver_options to be the current eigenvectors
        if (this->solver_options.get_solver_type() == numopt::eigenproblem::SolverType::DAVIDSON) {
            auto davidson_solver_options = dynamic_cast<numopt::eigenproblem::DavidsonSolverOptions&>(this->solver_options);

            for (size_t i = 0; i < this->solver_options.number_of_requested_eigenpairs; i++) {
                davidson_solver_options.X_0.col(i) = this->doci_solver.get_wavefunction(i).get_coefficients();
            }

            this->solver_options = davidson_solver_options;
        }
    }  // while not converged
}


/**
 *  @return a WaveFunction instance after performing the orbital optimization for a given eigenvector at @param index
 */
GQCG::WaveFunction DOCINewtonOrbitalOptimizer::get_wavefunction(size_t index) {
    if (this->is_converged) {
        return this->doci_solver.get_wavefunction(index);
    } else {
        throw std::logic_error("You are trying to get eigenpairs but the orbital optimization hasn't converged (yet).");
    }
}


}  // namespace GQCG
