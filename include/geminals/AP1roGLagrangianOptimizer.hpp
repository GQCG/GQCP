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
#ifndef AP1roGLagrangianOptimizer_hpp
#define AP1roGLagrangianOptimizer_hpp


#include "Geminals/BaseAP1roGSolver.hpp"
#include "Geminals/AP1roGVariables.hpp"


namespace GQCP {


/**
 *  A class that is able to optimize the AP1roG PSE Lagrangian
 */
class AP1roGLagrangianOptimizer : public BaseAP1roGSolver {
private:
    AP1roGVariables multipliers;  // the Lagrangian multipliers

public:
    // CONSTRUCTORS
    using BaseAP1roGSolver::BaseAP1roGSolver;  // inherit base constructors


    // GETTERS
    const AP1roGVariables& get_multipliers() const { return this->multipliers; }


    // PUBLIC METHODS
    void solve() override;
};


}  // namespace GQCP


#endif  /* AP1roGLagrangianOptimizer_hpp */
