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


#include "HamiltonianParameters/AtomicDecompositionParameters.hpp"

namespace GQCP {



/**
 *  @param molecular_hamiltonian_parameters     the complete molecular Hamiltonian parameters
 *  @param net_atomic_parameters                collection of net atomic Hamiltonian parameters
 *  @param interaction_parameters               collection of atomic interaction Hamiltonian parameters
 *  @param atomic_parameters                    collection of atomic Hamiltonian parameters
 */
AtomicDecompositionParameters::AtomicDecompositionParameters (const HamiltonianParameters<double>& molecular_hamiltonian_parameters, const std::vector<HamiltonianParameters<double>>& net_atomic_parameters,
                               const std::vector<HamiltonianParameters<double>>& interaction_parameters, const std::vector<HamiltonianParameters<double>>& atomic_parameters) :
        molecular_hamiltonian_parameters (molecular_hamiltonian_parameters),
        net_atomic_parameters (net_atomic_parameters),
        interaction_parameters (interaction_parameters),
        atomic_parameters (atomic_parameters)
{}


/**
 *  Constructs net atomic, atomic and atomic interaction Hamiltonian parameters in the AO basis for a diatomic molecule AB.
 *   the term "Nuclear" concerns how the electronic nuclear integrals (potential energy) are decomposed. The potential energy
 *   for basis functions on atom A for the charge on B are included in the interaction energy and not in the net atomic energy.
 *
 *  @param molecule     the molecule for which the AtomicDecompositionParameters should be calculated
 *  @param basisset     the name of the basisset corresponding to the AO basis
 *
 *  @return Atomic decomposed parameters:
 *      - net atomic parameters, HamiltonianParameters with:
 *          - one-electron nuclear integrals separated by atomic core and the atomic basis functions centered on that atom.
 *          - one-electron kinetic integrals separated per set of atomic basis functions centered on an atom.
 *          - two-electron integrals separated per set of atomic basis functions centered on an atom.
 *      - interaction parameters, HamiltonianParameters with:
 *          - remaining one- and two-electron contributions when deducting the net atomic parameters from the total HamiltonianParameters
 *          - scalar : nuclear repulsion
 *      - atomic parameters, HamiltonianParameters with:
 *          - net atomic parameters + interaction parameters/2
 *
 *  Ordering of the atomic Hamiltonian parameters are dependant on the ordering of the atoms in the molecule
 *   for the molecule AB:
 *      net_atomic_parameters will contains parameters for A then B.
 *      interaction_parameters will contain parameters for the AB interaction.
 *      atomic_parameters will contain parameters for A then B.
 */
AtomicDecompositionParameters AtomicDecompositionParameters::Nuclear(const Molecule& molecule, const std::string& basisset_name) {

    auto atoms = molecule.get_atoms();
    auto ao_basis = std::make_shared<AOBasis>(molecule, basisset_name);

    auto K = ao_basis->numberOfBasisFunctions();
    SquareMatrix<double> T_total = SquareMatrix<double>::Identity(K, K);

    if (atoms.size() > 2) {
        throw std::invalid_argument("AtomicDecompositionParameters::Nuclear(Molecule, std::string): Only available for diatomic molecules");
    }

    // retrieve the AObasis for the individual atoms so that we can retrieve net atomic nuclear integrals.
    AOBasis ao_basis_a ({{atoms[0]}}, basisset_name);
    AOBasis ao_basis_b ({{atoms[1]}}, basisset_name);

    OneElectronOperator<double> V_a = ao_basis_a.calculateLibintNuclearIntegrals();
    OneElectronOperator<double> V_b = ao_basis_b.calculateLibintNuclearIntegrals();

    // T_a and T_b are equal to the corresponding block from the molecular kinetic integrals (T_a = T.block(0,0, K_a, K_a))
    OneElectronOperator<double> T_a = ao_basis_a.calculateLibintKineticIntegrals();
    OneElectronOperator<double> T_b = ao_basis_b.calculateLibintKineticIntegrals();

    auto K_a = ao_basis_a.numberOfBasisFunctions();
    auto K_b = ao_basis_b.numberOfBasisFunctions();

    // create partition matrices for both atoms
    auto p_a = SquareMatrix<double>::PartitionMatrix(0, K_a, K);
    auto p_b = SquareMatrix<double>::PartitionMatrix(K_a, K_b, K);

    // retrieve the molecular integrals
    const auto& S = ao_basis->calculateLibintOverlapIntegrals();
    const auto& T = ao_basis->calculateLibintKineticIntegrals();
    const auto& V = ao_basis->calculateLibintNuclearIntegrals();
    const auto& g = ao_basis->calculateLibintCoulombRepulsionIntegrals();
    auto repulsion = molecule.calculateInternuclearRepulsionEnergy();

    OneElectronOperator<double> H = T + V;

    // Decompose the integrals corresponding to the formula's in Mario's thesis
    OneElectronOperator<double> h_a = OneElectronOperator<double>::Zero(K, K);
    OneElectronOperator<double> h_b = OneElectronOperator<double>::Zero(K, K);

    h_a.block(0, 0, K_a, K_a) = T_a + V_a;
    h_b.block(K_a , K_a, K_b, K_b) = T_b + V_b;

    OneElectronOperator<double> h_ab = H - h_a - h_b;

    auto g_a = g;
    auto g_b = g;
    auto g_ab = g;
    auto g_ba = g;

    g_a.matrixContraction<double>(p_a, 0);
    g_a.matrixContraction<double>(p_a, 2);

    g_b.matrixContraction<double>(p_b, 0);
    g_b.matrixContraction<double>(p_b, 2);

    g_ab.matrixContraction<double>(p_a, 0);
    g_ab.matrixContraction<double>(p_b, 2);

    g_ba.matrixContraction<double>(p_b, 0);
    g_ba.matrixContraction<double>(p_a, 2);

    GQCP::TwoElectronOperator<double> g_abba = g_ab.Eigen() + g_ba.Eigen();

    HamiltonianParameters<double> HAA(ao_basis, S, h_a, g_a, T_total);
    HamiltonianParameters<double> HBB(ao_basis, S, h_b, g_b, T_total);
    HamiltonianParameters<double> HAB(ao_basis, S, h_ab, g_abba, T_total, repulsion);
    HamiltonianParameters<double> HA(ao_basis, S, h_a + h_ab/2, g_a.Eigen() + (0.5)*g_abba.Eigen(), T_total, repulsion/2);
    HamiltonianParameters<double> HB(ao_basis, S, h_b + h_ab/2, g_b.Eigen() + (0.5)*g_abba.Eigen(), T_total, repulsion/2);

    std::vector<HamiltonianParameters<double>> net_atomic_parameters = {HAA, HBB};
    std::vector<HamiltonianParameters<double>> interaction_parameters = {HAB};
    std::vector<HamiltonianParameters<double>> atomic_parameters = {HA, HB};

    return AtomicDecompositionParameters(HamiltonianParameters<double>::Molecular(molecule, basisset_name), net_atomic_parameters, interaction_parameters, atomic_parameters);
}


}  // namespace GQCP
