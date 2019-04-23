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
#ifndef GQCP_ATOMICDECOMPOSITIONPARAMETERS_HPP
#define GQCP_ATOMICDECOMPOSITIONPARAMETERS_HPP



namespace GQCP {


/**
 *  A struct that holds the atomic decomposed Hamiltonian Parameters.
 *
 *  @tparam Scalar      the scalar type
 */
template <typename Scalar>
struct AtomicDecompositionParameters {

    Molecule molecule;  // decomposed molecule

    std::vector<HamiltonianParameters<double>> net_atomic_parameters;  // vector of net atomic Hamiltonian parameters, E_AA, E_BB
    std::vector<HamiltonianParameters<double>> interaction_parameters;  // vector of interaction Hamiltonian parameters, E_AB
    std::vector<HamiltonianParameters<double>> fragment_parameters;  // vector of total atomic or fragment contributions, E_AA + E_AB/2, E_BB + E_AB/2

}  // namespace GQCP


#endif //GQCP_ATOMICDECOMPOSITIONPARAMETERS_HPP
