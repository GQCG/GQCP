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
#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "CISolver/CISolver.hpp"
#include "HamiltonianBuilder/Hubbard.hpp"
#include "RDM/RDMCalculator.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>


namespace GQCP {
/**
*  A class that is able to drive Hubbard calculations.
*/
class HubbardDriver {
private:
    std::shared_ptr<GQCP::ProductFockSpace> fock_space;
    std::shared_ptr<GQCP::Hubbard> hubbard;
    std::shared_ptr<GQCP::CISolver> solver;
    GQCP::DenseSolverOptions dense_solver_options;

public:
    HubbardDriver(std::string csline, size_t num_states, size_t num_alpha, size_t num_beta, size_t K);

    std::vector<double> get_energies();

    std::vector<Eigen::MatrixXd> get_first_order_rdms();
};
}
