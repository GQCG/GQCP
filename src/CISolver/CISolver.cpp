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
#include "CISolver/CISolver.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  Constructor given a @param hamiltonian_builder and @param hamiltonian_parameters
 */
CISolver::CISolver(HamiltonianBuilder& hamiltonian_builder, const HamiltonianParameters& hamiltonian_parameters) :
    hamiltonian_builder (&hamiltonian_builder),
    hamiltonian_parameters (hamiltonian_parameters)
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

            numopt::eigenproblem::DenseSolver solver = numopt::eigenproblem::DenseSolver(matrix, dynamic_cast<numopt::eigenproblem::DenseSolverOptions&>(solver_options));

            solver.solve();
            this->eigenpairs = solver.get_eigenpairs();

            break;
        }

        case numopt::eigenproblem::SolverType::DAVIDSON: {

            Eigen::VectorXd diagonal = this->hamiltonian_builder->calculateDiagonal(this->hamiltonian_parameters);
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
}

/**
 *  @return WaveFunction instance after solving the CI problem for a given eigenvector at @param index
 */
GQCG::WaveFunction CISolver::get_wavefunction(size_t index) {
    if (index > this->eigenpairs.size()) {
        throw std::logic_error("Not enough requested eigenpairs for the given index.");
    }
    return WaveFunction(*this->hamiltonian_builder->get_fock_space(), this->eigenpairs[index].get_eigenvector());
}


}  // namespace GQCP
