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


#include "Mathematical/Optimization/Eigenproblem/EigenproblemEnvironment.hpp"
#include "ONVBasis/SpinResolvedFrozenONVBasis.hpp"
#include "ONVBasis/SpinResolvedONVBasis.hpp"
#include "ONVBasis/SpinResolvedSelectedONVBasis.hpp"
#include "ONVBasis/SpinUnresolvedFrozenONVBasis.hpp"
#include "ONVBasis/SpinUnresolvedONVBasis.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"


namespace GQCP {
namespace CIEnvironment {


/**
 *  @return an environment suitable for solving spin-resolved FCI eigenvalue problems.
 */
template <typename Scalar>
EigenproblemEnvironment Dense(const SQHamiltonian<Scalar>& sq_hamiltonian, const SpinResolvedONVBasis& onv_basis) {

    const auto H = onv_basis.evaluateOperatorDense(sq_hamiltonian, true);  // 'true': with calculation of the diagonal
    return EigenproblemEnvironment::Dense(H);
}


}  // namespace CIEnvironment
}  // namespace GQCP
