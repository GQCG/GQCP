// This file is part of GQCG-GQCP.
//
// Copyright (C) 2017-2020  the GQCG developers
//
// GQCG-GQCP is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GQCG-GQCP is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-GQCP.  If not, see <http://www.gnu.org/licenses/>.

#pragma once


#include "Mathematical/Algorithm/Algorithm.hpp"
#include "Mathematical/Algorithm/IterativeAlgorithm.hpp"
#include "Mathematical/Optimization/Eigenproblem/EigenproblemEnvironment.hpp"
#include "ONVBasis/BaseONVBasis.hpp"
#include "QCMethod/QCStructure.hpp"
#include "QCModel/CI/LinearExpansion.hpp"

#include <memory>


namespace GQCP {
namespace QCMethod {


/**
 *  The configuration interaction quantum chemical method.
 * 
 *  @tparam _ONVBasis           the type of ONV basis
 */
template <typename _ONVBasis>
class CI {

public:
    using ONVBasis = _ONVBasis;


private:
    size_t number_of_states;  // the number of states that searched for (incuding the ground state)
    ONVBasis onv_basis;


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param onv_basis                        the ONV basis with respect to which the configuration interaction is expressed
     *  @param number_of_states                 the number of states that searched for (including the ground state)
     */
    CI(const ONVBasis& onv_basis, const size_t number_of_states = 1) :
        onv_basis {onv_basis},
        number_of_states {number_of_states} {}


    /*
     *  PUBLIC METHODS
     */

    /**
     *  Optimize the electronic structure model.
     * 
     *  @tparam Solver              the type of the solver
     * 
     *  @param solver               the solver that will try to optimize the parameters
     */
    template <typename Solver>
    QCStructure<LinearExpansion<ONVBasis>> optimize(Solver& solver, EigenproblemEnvironment& environment) const {

        // The CI method's responsibility is to try to optimize the parameters of its method, given a solver and associated environment.
        solver.perform(environment);


        // Extract the requested number of eigenpairs from the environment and place them into the LinearExpansion wave function model. Consequently, check if the LinearExpansion fulfills the objective.
        const auto eigenpairs = environment.eigenpairs(this->number_of_states);

        std::vector<LinearExpansion<ONVBasis>> linear_expansions {};
        linear_expansions.reserve(number_of_states);

        std::vector<double> energies {};
        energies.reserve(this->number_of_states);

        for (const auto& eigenpair : eigenpairs) {
            linear_expansions.emplace_back(onv_basis, eigenpair.eigenvector());

            // TODO: We can't check for a 'NormedHamiltonianEigenvectorObjective' to be fulfilled, because there is no uniform way to let an SQHamiltonian act on a LinearExpansion.

            energies.push_back(eigenpair.eigenvalue());
        }


        // Wrap all the requested number of states into a QCStructure.
        // Since we have already created a list of LinearExpansions, we only have to create a list of the corresponding energies.
        return QCStructure<LinearExpansion<ONVBasis>>(energies, linear_expansions);
    }
};


}  // namespace QCMethod
}  // namespace GQCP
