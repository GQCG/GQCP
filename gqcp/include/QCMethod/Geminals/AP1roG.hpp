// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2019  the GQCG developers
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
#pragma once


#include "Mathematical/Optimization/NonLinearEquation/NonLinearEquationEnvironment.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/QCStructure.hpp"
#include "QCModel/Geminals/AP1roG.hpp"


namespace GQCP {
namespace QCMethod {


/**
 *  The AP1roG quantum chemical method.
 */
class AP1roG {


private:
    size_t N_P;  // the number of electron pairs
    size_t K;  // the number of spatial orbitals

    SQHamiltonian<double> sq_hamiltonian;  // the second-quantized Hamiltonian in an orthonormal basis


public:

    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param sq_hamiltonian           the second-quantized Hamiltonian in an orthonormal basis
     *  @param N_P                      the number of electron pairs
     */
    AP1roG(const SQHamiltonian<double>& sq_hamiltonian, const size_t N_P) :
        N_P (N_P),
        K (sq_hamiltonian.dimension()),  // number of spatial orbitals
        sq_hamiltonian (sq_hamiltonian)
    {}


    /*
     *  PUBLIC METHODS
     */

    /**
     *  Optimize the AP1roG wave function model: find the parameters satisfy the given objective.
     * 
     *  @tparam Solver              the type of the solver
     * 
     *  @param solver               the solver that will try to optimize the parameters
     *  @param environment          the environment, which acts as a sort of calculation space for the solver
     */
    template <typename Solver>
    QCStructure<GQCP::QCModel::AP1roG> optimize(Solver& solver, NonLinearEquationEnvironment<double>& environment) const {

        // The AP1roG method's responsibility is to try to optimize the parameters of its method, given a solver and associated environment.
        solver.perform(environment);

        // To make a QCStructure, we need the electronic energy, geminal coefficients and number of electrons.
        // Furthermore, the solvers only find the ground state wave function parameters, so the QCStructure only needs to contain the parameters for one state.
        const auto G_optimal = AP1roGGeminalCoefficients::FromColumnMajor(environment.variables.back(), this->N_P, this->K);
        const GQCP::QCModel::AP1roG ap1rog_parameters {G_optimal};

        const auto E_electronic = ap1rog_parameters.calculateEnergy(sq_hamiltonian);

        return QCStructure<GQCP::QCModel::AP1roG>({E_electronic}, {ap1rog_parameters});
    }
};


}  // namespace QCMethod
}  // namespace GQCP
