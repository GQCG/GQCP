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
#ifndef GQCP_QCMETHOD_HUBBARD_HPP
#define GQCP_QCMETHOD_HUBBARD_HPP


#include "RDM/OneRDM.hpp"


namespace GQCP {
namespace QCMethod {


/**
 *  A class that is a wrapper around solving the dense eigenvalue problem for Hubbard
 */
class Hubbard {
private:
    std::string csline;  // a comma-separated line that contains the upper (or lower) triangle of the Hubbard hopping matrix
    size_t num_states;  // the number of states that should be targeted (the lowest num_states are found)
    size_t N_alpha;  // the number of alpha electrons
    size_t N_beta;  // the number of beta electrons
    size_t K;  // the number of spatial orbitals

    bool is_solved = false;
    std::vector<double> energy_solutions;  // the found energies
    std::vector<OneRDM<double>> dm_solutions;  // the found DMs


public:
    // CONSTRUCTORS

    /**
     *  @param cslin            a comma-separated line that contains the upper (or lower) triangle of the Hubbard hopping matrix
     *  @param num_states       the number of states that should be targeted (the lowest num_states are found)
     *  @param num_alpha        the number of alpha electrons
     *  @param num_beta         the number of beta electrons
     *  @param num_orb          the number of spatial orbitals
     */
    Hubbard(const std::string& csline, const size_t num_states, const size_t num_alpha, const size_t num_beta);


    // PUBLIC METHODS

    /**
     *  Solve the Hubbard eigenvalue problem
     */
    void solve();

    /**
     *  @return the lowest (requested) energies
     */
    const std::vector<double>& energies() const;

    /**
     *  @return the DMs that correspond to the lowest (requested) wave functions
     */
    const std::vector<OneRDM<double>>& oneRDMs() const;
};


}  // namespace QCMethod
}  // namespace GQCP


#endif  // GQCP_QCMETHOD_HUBBARD_HPP
