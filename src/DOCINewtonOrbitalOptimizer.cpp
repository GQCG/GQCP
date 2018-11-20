// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2018  the GQCG developers
// 
// GQCG-gqcp is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// GQCG-gqcp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-gqcp.  If not, see <http://www.gnu.org/licenses/>.
// 
#include "DOCINewtonOrbitalOptimizer.hpp"

#include <cpputil.hpp>
#include <unsupported/Eigen/MatrixFunctions>

#include "CISolver/CISolver.hpp"
#include "RDM/RDMCalculator.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */
/**
 *  @param doci         the DOCI HamiltonianBuilder
 *  @param ham_par      the Hamiltonian parameters in an orthonormal basis
 */
DOCINewtonOrbitalOptimizer::DOCINewtonOrbitalOptimizer(const GQCP::DOCI& doci, const GQCP::HamiltonianParameters& ham_par) :
    doci (doci),
    ham_par (ham_par)
{}


/*
 *  GETTERS
 */
const std::vector<numopt::eigenproblem::Eigenpair>& DOCINewtonOrbitalOptimizer::get_eigenpairs() const {
    if (this->is_converged) {
        return this->eigenpairs;
    } else {
        throw std::logic_error("You are trying to get eigenpairs but the orbital optimization hasn't converged (yet).");
    }
}

const numopt::eigenproblem::Eigenpair& DOCINewtonOrbitalOptimizer::get_eigenpair(size_t index) const {
    if (this->is_converged) {
        return this->eigenpairs[index];
    } else {
        throw std::logic_error("You are trying to get eigenpairs but the orbital optimization hasn't converged (yet).");
    }
}



/*
 *  PUBLIC METHODS
 */
/**
 *  Do the orbital optimization for DOCI
 *
 *  @param solver_options       solver options for the CI solver
 *  @param oo_options           options for the orbital optimization
 */
void DOCINewtonOrbitalOptimizer::solve(numopt::eigenproblem::BaseSolverOptions& solver_options, const GQCP::OrbitalOptimizationOptions& oo_options) {
    this->is_converged = false;
    auto K = this->ham_par.get_K();

    size_t oo_iterations = 0;
    while (!(this->is_converged)) {

        // Solve the DOCI eigenvalue equation, using the options provided
        GQCP::CISolver doci_solver (this->doci, this->ham_par);  // update the CI solver with the rotated Hamiltonian parameters
        doci_solver.solve(solver_options);
        GQCP::WaveFunction ground_state = doci_solver.get_wavefunction();

        // Calculate the 1- and 2-RDMs
        GQCP::RDMCalculator rdm_calculator (*(this->doci.get_fock_space()));
        auto D = rdm_calculator.calculate1RDMs(ground_state.get_coefficients()).one_rdm;  // spin-summed 1-RDM
        auto d = rdm_calculator.calculate2RDMs(ground_state.get_coefficients()).two_rdm;  // spin-summed 2-RDM


        // Calculate the electronic gradient at kappa = 0
        Eigen::MatrixXd F = this->ham_par.calculateGeneralizedFockMatrix(D, d).get_matrix_representation();
        Eigen::MatrixXd gradient_matrix = 2 * (F - F.transpose());
        Eigen::VectorXd gradient_vector = cpputil::linalg::strictLowerTriangle(gradient_matrix);  // gradient vector with the free parameters, at kappa = 0


        // Calculate the electronic Hessian at kappa = 0
        GQCP::TwoElectronOperator W = this->ham_par.calculateSuperGeneralizedFockMatrix(D, d);
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
        numopt::VectorFunction gradient_function = [gradient_vector](const Eigen::VectorXd& x) { return gradient_vector; };
        numopt::MatrixFunction hessian_function = [hessian_matrix](const Eigen::VectorXd& x) { return hessian_matrix; };

        Eigen::VectorXd kappa_vector = numopt::newtonStep(Eigen::VectorXd::Zero(K), gradient_function, hessian_function);  // with only the free parameters


        // If the calculated norm is zero, we have reached a critical point
        if (gradient_vector.norm() < oo_options.convergence_threshold) {

            // If we have found a critical point, but we have a negative eigenvalue for the Hessian, continue in that direction
            if (hessian_solver.eigenvalues()(0) < 0) {
                kappa_vector = hessian_solver.eigenvectors().col(0);
            }
            else {  // the Hessian is confirmed to be positive definite, so we have reached a minimum
                this->is_converged = true;

                // Set solutions
                this->eigenpairs = doci_solver.get_eigenpairs();

                break;  // no need to continue if we have converged
            }


        } else {
            oo_iterations++;

            if (oo_iterations >= oo_options.maximum_number_of_iterations) {
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
        if (solver_options.get_solver_type() == numopt::eigenproblem::SolverType::DAVIDSON) {
            auto davidson_solver_options = dynamic_cast<numopt::eigenproblem::DavidsonSolverOptions&>(solver_options);

            for (size_t i = 0; i < solver_options.number_of_requested_eigenpairs; i++) {
                davidson_solver_options.X_0.col(i) = doci_solver.get_wavefunction(i).get_coefficients();
            }

            solver_options = davidson_solver_options;
        }
    }  // while not converged
}


/**
 *  @param index        the index of the index-th excited state
 *
 *  @return the index-th excited state after doing the OO-DOCI calculation
 */
GQCP::WaveFunction DOCINewtonOrbitalOptimizer::get_wavefunction(size_t index) {
    if (index > this->eigenpairs.size()) {
        throw std::logic_error("Not enough requested eigenpairs for the given index.");
    }
    return WaveFunction(*dynamic_cast<GQCP::FockSpace*>(this->doci.get_fock_space()), this->eigenpairs[index].get_eigenvector());
}


}  // namespace GQCP
