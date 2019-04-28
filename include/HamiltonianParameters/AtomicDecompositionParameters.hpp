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


#include <Molecule.hpp>
#include "HamiltonianParameters/HamiltonianParameters.hpp"

namespace GQCP {


/**
 *  A struct that holds the collection of Hamiltonian Parameters that represent different molecular decompositions.
 *
 *  Currently only implemented for diatomic molecules
 */
struct AtomicDecompositionParameters {

    Molecule molecule;  // Decomposed molecule

    HamiltonianParameters<double> molecular_hamiltonian_parameters;  // the Hamiltonian Parameters of the complete molecule
    std::vector<HamiltonianParameters<double>> net_atomic_parameters;  // vector of net atomic Hamiltonian parameters, E_AA, E_BB
    std::vector<HamiltonianParameters<double>> interaction_parameters;  // vector of interaction Hamiltonian parameters, E_AB
    std::vector<HamiltonianParameters<double>> fragment_parameters;  // vector of total atomic or fragment contributions, E_AA + E_AB/2, E_BB + E_AB/2


    explicit AtomicDecompositionParameters (const Molecule& molecule, const std::string& basisset_name) : molecule(molecule), molecular_hamiltonian_parameters(HamiltonianParameters<double>::Molecular(molecule, basisset_name)) {

        auto atoms = molecule.get_atoms();

        if (atoms.size() > 2) {
            throw std::invalid_argument("HamiltonianParameters::atomicDecomposition(): The Hamiltonian parameters are set up for more than 2 atoms, currently only available for diatomic molecules");
        }

        Molecule atom_a ({atoms[0]});
        Molecule atom_b ({atoms[1]});

        auto ao_basis = std::make_shared<AOBasis>(molecule, basisset_name);
        AOBasis ao_basis_a (atom_a, basisset_name);
        AOBasis ao_basis_b (atom_b, basisset_name);

        auto Ka = ao_basis_a.numberOfBasisFunctions();
        auto Kb = ao_basis_b.numberOfBasisFunctions();
        auto K = ao_basis->numberOfBasisFunctions();

        auto p_a = SquareMatrix<double>::PartitionMatrix(0, Ka, K);
        auto p_b = SquareMatrix<double>::PartitionMatrix(Ka, Kb, K);

        SquareMatrix<double> identity = SquareMatrix<double>::Identity(K, K);

        const auto& S = ao_basis->calculateOverlapIntegrals();
        const auto& T = ao_basis->calculateKineticIntegrals();
        const auto& V = ao_basis->calculateNuclearIntegrals();
        const auto& g = ao_basis->calculateCoulombRepulsionIntegrals();
        auto repulsion = molecule.calculateInternuclearRepulsionEnergy();

        OneElectronOperator<double> H = T + V;

        OneElectronOperator<double> h_a = p_a * H * p_a;
        OneElectronOperator<double> h_b = p_b * H * p_b;
        OneElectronOperator<double> h_ab = p_b * H * p_a + p_a * H * p_b;

        auto g_a = g;
        auto g_b = g;
        auto g_ab = g;
        auto g_ba = g;

        g_a.fourModeMultiplication<double>(p_a, identity, p_a, identity);
        g_b.fourModeMultiplication<double>(p_b, identity, p_b, identity);
        g_ab.fourModeMultiplication<double>(p_a, identity, p_b, identity);
        g_ba.fourModeMultiplication<double>(p_b, identity, p_a, identity);
        GQCP::TwoElectronOperator<double> g_abba = g_ab.Eigen() + g_ba.Eigen();

        HamiltonianParameters<double> HAA(ao_basis, S, h_a, g_a, identity);
        HamiltonianParameters<double> HBB(ao_basis, S, h_b, g_b, identity);
        HamiltonianParameters<double> HAB(ao_basis, S, h_ab, g_abba, identity, repulsion);
        HamiltonianParameters<double> HA(ao_basis, S, h_a + h_ab/2, g_a.Eigen() + (0.5)*g_abba.Eigen(), identity);
        HamiltonianParameters<double> HB(ao_basis, S, h_b + h_ab/2, g_b.Eigen() + (0.5)*g_abba.Eigen(), identity);

        this->net_atomic_parameters = {HAA, HBB};
        this->interaction_parameters = {HAB};
        this->fragment_parameters = {HA, HB};
    };
};

}  // namespace GQCP


#endif //GQCP_ATOMICDECOMPOSITIONPARAMETERS_HPP
