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


#include "Mathematical/Optimization/Eigenproblem/GeneralizedEigenproblemEnvironment.hpp"
#include "Molecule/Molecule.hpp"
#include "Operator/FirstQuantized/NuclearRepulsionOperator.hpp"


namespace GQCP {
namespace NOCIEnvironment {


/**
 *  Create an environment suitable for solving dense NOCI eigenvalue problems for the given Hamiltonian and non-orthogonal basis.
 *
 *  @tparam Hamiltonian                       The type of Hamiltonian whose eigenproblem is trying to be solved.
 *  @tparam NonOrthogonalBasis                The type of non-orthogonal basis in which the Hamiltonian should be represented.
 *
 *  @param hamiltonian                        A second-quantized Hamiltonian expressed in an orthonormal orbital basis.
 *  @param non_orthogonal_basis               A non-orthogonal basis in which the Hamiltonian eigenproblem should be solved.
 *  @param molecule                           The molecule for which the eigenproblem will be solved. The NOCI Hamiltonian requires a nuclear repulsion component, which can be calculated from the molecule.
 *
 *  @return An `EigenproblemEnvironment` initialized suitable for solving dense NOCI eigenvalue problems for the given Hamiltonian and non-orthogonal basis.
 */
template <typename Hamiltonian, typename NonOrthogonalBasis>
auto Dense(const Hamiltonian& hamiltonian, const NonOrthogonalBasis& non_orthogonal_basis, const Molecule& molecule) -> GeneralizedEigenproblemEnvironment<typename Hamiltonian::Scalar> {

    using Scalar = typename Hamiltonian::Scalar;  // The scalar type of a Hamiltonian element.

    // Determine the dense matrix representation of the Hamiltonian in the given non-orthogonal basis, and supply it to an `EigenproblemEnvironment`.
    const auto H = non_orthogonal_basis.evaluateHamiltonianOperator(hamiltonian, NuclearRepulsionOperator(molecule.nuclearFramework()));
    const auto S = non_orthogonal_basis.evaluateOverlapOperator();
    return GeneralizedEigenproblemEnvironment<Scalar>::Dense(H, S);
}


}  // namespace NOCIEnvironment
}  // namespace GQCP
