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


#include "Mathematical/Algorithm/Step.hpp"
#include "Mathematical/Optimization/Eigenproblem/EigenproblemEnvironment.hpp"


namespace GQCP {


class GuessVectorUpdate :
    public Step<EigenproblemEnvironment> {

public:

    /*
     *  OVERRIDDEN PUBLIC METHODS
     */

    /**
     *  Calculate new guesses for the eigenvectors from the diagonalized subspace matrix
     * 
     *  @param environment              the environment that acts as a sort of calculation space
     */
    void execute(EigenproblemEnvironment& environment) override {

        // X contains the new guesses for the eigenvectors, V is the subspace and Z are the eigenvectors of the subspace matrix
        environment.X = environment.V * environment.Z;  // X is a linear combination of the current subspace vectors
        environment.eigenvectors = environment.X;
    }
};


}  // namespace GQCP
