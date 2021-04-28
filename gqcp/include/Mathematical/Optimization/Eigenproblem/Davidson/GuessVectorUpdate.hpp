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


#include "Mathematical/Algorithm/Step.hpp"
#include "Mathematical/Optimization/Eigenproblem/EigenproblemEnvironment.hpp"


namespace GQCP {


/**
 *  A step that calculates new guesses for the eigenvectors from the diagonalized subspace matrix.
 */
class GuessVectorUpdate:
    public Step<EigenproblemEnvironment<double>> {

public:
    /*
     *  PUBLIC OVERRIDDEN METHODS
     */

    /**
     *  @return a textual description of this algorithmic step
     */
    std::string description() const override {
        return "Calculate new guesses for the eigenvectors from the diagonalized subspace matrix.";
    }


    /**
     *  Calculate new guesses for the eigenvectors from the diagonalized subspace matrix.
     * 
     *  @param environment              the environment that acts as a sort of calculation space
     */
    void execute(EigenproblemEnvironment<double>& environment) override {

        // X contains the new guesses for the eigenvectors, V is the subspace and Z are the eigenvectors of the subspace matrix.
        environment.X = environment.V * environment.Z;  // X is a linear combination of the current subspace vectors
        environment.eigenvectors = environment.X;
    }
};


}  // namespace GQCP
