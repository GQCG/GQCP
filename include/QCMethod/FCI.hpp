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
#include "HamiltonianBuilder/FCI.hpp"
#include "FockSpace/ProductFockSpace.hpp"
#include "RDM/RDMCalculator.hpp"


namespace GQCP {
/**
*  A class that is able to drive FCI calculations.
*/
class FullConfigurationInteractionDriver {
private:
    std::shared_ptr<GQCP::Molecule> molecule;
    std::shared_ptr<GQCP::CISolver> solver;

public:
    FullConfigurationInteractionDriver(std::string xyz_filename, std::string basis_set, size_t num_alpha, size_t num_beta);

    double get_energy();
};
}
