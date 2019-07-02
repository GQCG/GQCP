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
#include "QCMethod/Hubbard.hpp"

#include "CISolver/CISolver.hpp"
#include "HamiltonianBuilder/Hubbard.hpp"
#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "RDM/RDMCalculator.hpp"

#include <boost/algorithm/string.hpp>

#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>


namespace GQCP {
namespace QCMethod {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param cslin            a comma-separated line that contains the upper (or lower) triangle of the Hubbard hopping matrix
 *  @param num_states       the number of states that should be targeted (the lowest num_states are found)
 *  @param num_alpha        the number of alpha electrons
 *  @param num_beta         the number of beta electrons
 *  @param num_orb          the number of spatial orbitals
 */
Hubbard::Hubbard(const std::string& csline, const size_t num_states, const size_t num_alpha, const size_t num_beta, const size_t num_orb) :
    csline (csline),
    num_states (num_states),
    N_alpha (num_alpha),
    N_beta (num_beta),
    K (num_orb)
{}



/*
 *  PUBLIC METHODS
 */

/**
 *  Solve the Hubbard eigenvalue problem
 */
void Hubbard::solve() {

    std::string triagonal_line = csline;

    // Read the upper triagonal of the hopping matrix
    GQCP::VectorX<double> triagonal;
    if (csline.empty()) {
        throw std::invalid_argument("hubbard(driver): Comma-separated was empty!");
    }

    // Split comma-separated line
    std::vector<std::string> splitted_line;
    boost::split(splitted_line, triagonal_line, boost::is_any_of(","));

    std::vector<double> triagonal_data;
    for (const std::string &x : splitted_line) {
        triagonal_data.push_back(std::stod(x));
    }

    triagonal = Eigen::Map<Eigen::VectorXd>(triagonal_data.data(), triagonal_data.size());
    GQCP::HoppingMatrix H = GQCP::HoppingMatrix::FromUpperTriangle(triagonal);


    // Actual calculations
    const auto ham_par = GQCP::HamiltonianParameters<double>::Hubbard(H);
    if (ham_par.get_K() != K) {
        throw std::invalid_argument("hubbard(driver): The given number of sites does not match the triagonal");
    }


    // Initialize and solve the Hubbard eigenvalue problem
    ProductFockSpace fock_space (this->K, this->N_alpha, this->N_beta);
    GQCP::Hubbard hubbard_builder (fock_space);
    CISolver ci_solver (hubbard_builder, ham_par);

    DenseSolverOptions ci_solver_options;
    ci_solver_options.number_of_requested_eigenpairs = this->num_states;
    ci_solver.solve(ci_solver_options);
    this->is_solved = true;


    // Store the solutions
    GQCP::RDMCalculator rdm_calculator (fock_space);
    for (const auto& eigenpair : ci_solver.get_eigenpairs()) {
        this->energy_solutions.push_back(eigenpair.get_eigenvalue());
        
        rdm_calculator.set_coefficients(eigenpair.get_eigenvector());
        this->dm_solutions.push_back(rdm_calculator.calculate1RDMs().one_rdm);
    }
}


/**
 *  @return the lowest (requested) energies
 */
const std::vector<double>& Hubbard::energies() const {

    if (this->is_solved) {
        return this->energy_solutions;
    } else {
        throw std::runtime_error("Hubbard::energies(): You are trying to get energies but the method hasn't been solved yet.");
    }
}


/**
 *  @return the DMs that correspond to the lowest (requested) wave functions
 */
const std::vector<OneRDM<double>>& Hubbard::oneRDMs() const {

    if (this->is_solved) {
        return this->dm_solutions;
    } else {
        throw std::runtime_error("Hubbard::oneRDMs(): You are trying to get DMs but the method hasn't been solved yet.");
    }
}


}  // namespace QCMethod
}  // namespace GQCP
