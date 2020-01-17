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


#include "Basis/TransformationMatrix.hpp"
#include "Mathematical/Optimization/IterativeSolver.hpp"


namespace GQCP {


/**
 *  A restricted Hartree-Fock self-consistent field (RHF SCF) solver.
 * 
 *  @tparam _ExpansionScalar        the type of scalar that is used to describe the expansion coefficients
 */
template <typename _ExpansionScalar>
class RHFSCFSolver : public IterativeSolver<TransformationMatrix<_ExpansionScalar>, SimpleSCFSolver<_ExpansionScalar>> {

public:

    /*
     *  CONSTRUCTORS
     */

    RHFSCFSolver


};


}  // namespace GQCP
