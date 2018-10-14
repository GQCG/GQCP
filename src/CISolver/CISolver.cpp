#include "CISolver/CISolver.hpp"


namespace GQCG {


/*
 *  CONSTRUCTORS
 */

/**
 *  Constructor given a @param hamiltonian_builder and @param hamiltonian_parameters
 */
CISolver::CISolver(HamiltonianBuilder& hamiltonian_builder, HamiltonianParameters& hamiltonian_parameters) :
    hamiltonian_builder(&hamiltonian_builder),
    hamiltonian_parameters(hamiltonian_parameters)
{
    auto K = hamiltonian_parameters.get_h().get_dim();
    if (K != this->hamiltonian_builder->get_fock_space()->get_K()) {
        throw std::invalid_argument("Basis functions of the Fock space and hamiltonian_parameters are incompatible.");
    }
}



/*
 *  CONSTRUCTORS
 */

/**
 *  solves the CI problem, setting the eigenvectors and -values
 */
void CISolver::solve(numopt::eigenproblem::BaseSolverOptions& solver_options) {
    numopt::eigenproblem::SolverType solver_type = solver_options.get_solver_type();
    switch (solver_type) {

        case numopt::eigenproblem::SolverType::DENSE: {
            Eigen::MatrixXd matrix = this->hamiltonian_builder->constructHamiltonian(this->hamiltonian_parameters);
            numopt::eigenproblem::DenseSolver solver = numopt::eigenproblem::DenseSolver(matrix, dynamic_cast<numopt::eigenproblem::DenseSolverOptions&>(solver_options)));
            solver.solve();
            this->eigenpairs = solver.get_eigenpairs();
            break;
        }

        case numopt::eigenproblem::SolverType::DAVIDSON: {

            Eigen::VectorXd diagonal = this->hamiltonian_builder->calculateDiagonal(this->hamiltonian_parameters);

            // Davidson Solver requires us to specify the macvec

            numopt::VectorFunction matrixVectorProduct = [this, &diagonal](const Eigen::VectorXd& x) { return hamiltonian_builder->matrixVectorProduct(hamiltonian_parameters, x, diagonal); };
            numopt::eigenproblem::DavidsonSolver solver = numopt::eigenproblem::DavidsonSolver(matrixVectorProduct, diagonal, dynamic_cast<numopt::eigenproblem::DavidsonSolverOptions&>(solver_options));
            solver.solve();
            this->eigenpairs = solver.get_eigenpairs();
            break;

        }

        case numopt::eigenproblem::SolverType::SPARSE: {
            throw std::invalid_argument("Sparse not implemented");
            break;
        }

    }
    this->eigenpairs = solver->get_eigenpairs();
}

/**
 *  @return WaveFunction instance after solving the CI problem for a given eigenvector.
 */
WaveFunction CISolver::get_wavefunction(size_t index) {
    if (index < this->eigenpairs.size()) {
        throw std::logic_error("Non existent wave function, you requested to few eigenpairs in your solve (or did not solve)");
    }
    return WaveFunction(*this->hamiltonian_builder->get_fock_space(), this->eigenpairs[index].get_eigenvector());
}


}  // namespace GQCG