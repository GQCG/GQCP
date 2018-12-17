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

#include "optimization/DenseSolver.hpp"
#include "optimization/DavidsonSolver.hpp"
#include "optimization/SparseSolver.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param hamiltonian_builder      the HamiltonianBuilder for which the CI eigenvalue problem should be solved
 *  @param hamiltonian_parameters   the Hamiltonian parameters in an orthonormal basis
 */
CISolver::CISolver(const HamiltonianBuilder& hamiltonian_builder, const HamiltonianParameters& hamiltonian_parameters) :
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
 *  @param solver_options       specify a type of solver and its options
 *
 *  Solve the CI eigenvalue problem and set the eigenpairs internally
 */
void CISolver::solve(const BaseSolverOptions& solver_options) {

    switch (solver_options.get_solver_type()) {

        case SolverType::DENSE: {

            Eigen::MatrixXd matrix = this->hamiltonian_builder->constructHamiltonian(this->hamiltonian_parameters);

            DenseSolver solver (matrix, dynamic_cast<const DenseSolverOptions&>(solver_options));

            solver.solve();
            this->eigenpairs = solver.get_eigenpairs();

            break;
        }

        case SolverType::DAVIDSON: {

            Eigen::VectorXd diagonal = this->hamiltonian_builder->calculateDiagonal(this->hamiltonian_parameters);
            VectorFunction matrixVectorProduct = [this, &diagonal](const Eigen::VectorXd& x) { return hamiltonian_builder->matrixVectorProduct(hamiltonian_parameters, x, diagonal); };

            DavidsonSolver solver (matrixVectorProduct, diagonal, dynamic_cast<const DavidsonSolverOptions&>(solver_options));

            solver.solve();
            this->eigenpairs = solver.get_eigenpairs();

            break;
        }

        case SolverType::SPARSE: {
            throw std::invalid_argument("Sparse not implemented");
            break;
        }
    }
}


/**
 *  @param index        the index of the index-th excited state
 *
 *  @return the index-th excited state after solving the CI eigenvalue problem
 */
WaveFunction CISolver::makeWavefunction(size_t index) const {
    if (index > this->eigenpairs.size()) {
        throw std::logic_error("Not enough requested eigenpairs for the given index.");
    }
    return WaveFunction(*this->hamiltonian_builder->get_fock_space(), this->eigenpairs[index].get_eigenvector());
}


}  // namespace GQCP
