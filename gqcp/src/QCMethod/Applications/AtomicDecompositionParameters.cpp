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

#include "QCMethod/Applications/AtomicDecompositionParameters.hpp"

#include "Basis/SpinorBasis/RSpinorBasis.hpp"
#include "Operator/FirstQuantized/Operator.hpp"


namespace GQCP {


/**
 *  @param molecular_hamiltonian_parameters     the complete molecular Hamiltonian
 *  @param net_atomic_parameters                collection of the net atomic Hamiltonian
 *  @param interaction_parameters               collection of the atomic interaction Hamiltonian
 *  @param atomic_parameters                    collection of the atomic Hamiltonian
 */
AtomicDecompositionParameters::AtomicDecompositionParameters(const RSQHamiltonian<double>& molecular_hamiltonian_parameters, const std::vector<RSQHamiltonian<double>>& net_atomic_parameters, const std::vector<RSQHamiltonian<double>>& interaction_parameters, const std::vector<RSQHamiltonian<double>>& atomic_parameters) :
    molecular_hamiltonian {molecular_hamiltonian_parameters},
    net_atomic_parameters {net_atomic_parameters},
    interaction_parameters {interaction_parameters},
    atomic_parameters {atomic_parameters} {}


/**
 *  Constructs net atomic, atomic and atomic interaction Hamiltonian in the AO basis for a diatomic molecule AB.
 *   the term "Nuclear" concerns how the electronic nuclear integrals (potential energy) are decomposed. The potential energy
 *   for basis functions on atom A for the charge on B are included in the interaction energy and not in the net atomic energy.
 *
 *  @param molecule     the molecule for which the AtomicDecompositionParameters should be calculated
 *  @param basisset     the name of the basisset corresponding to the AO basis
 *
 *  @return Atomic decomposed parameters:
 *      - net atomic, with:
 *          - one-electron nuclear integrals separated by atomic core and the atomic basis functions centered on that atom.
 *          - one-electron kinetic integrals separated per set of atomic basis functions centered on an atom.
 *          - two-electron integrals separated per set of atomic basis functions centered on an atom.
 *      - interaction, with:
 *          - remaining one- and two-electron contributions when deducting the net atomic Hamiltonian from the total HamiltonianParameters
 *          - scalar : nuclear repulsion
 *      - atomic, with:
 *          - net atomic parameters + interaction parameters/2
 *
 *  Ordering of the atomic Hamiltonian are dependant on the ordering of the atoms in the molecule
 *   for the molecule AB:
 *      net_atomic_parameters will contain first A, then B.
 *      interaction_parameters will contain the AB interaction.
 *      atomic_parameters will contain first A, then B.
 */
AtomicDecompositionParameters AtomicDecompositionParameters::Nuclear(const Molecule& molecule, const std::string& basisset_name) {

    const auto& atoms = molecule.nuclearFramework().nucleiAsVector();
    if (atoms.size() > 2) {
        throw std::invalid_argument("AtomicDecompositionParameters::Nuclear(Molecule, std::string): Only available for diatomic molecules");
    }

    const RSpinorBasis<double, GTOShell> spinor_basis {molecule, basisset_name};
    const auto K = spinor_basis.numberOfSpatialOrbitals();


    // Retrieve an AO basis for the individual atoms so that we can retrieve net atomic nuclear integrals
    const NuclearFramework nuclear_framework_a {{atoms[0]}};
    const NuclearFramework nuclear_framework_b {{atoms[1]}};

    RSpinorBasis<double, GTOShell> spinor_basis_a {nuclear_framework_a, basisset_name};  // in non-orthogonal AO basis
    RSpinorBasis<double, GTOShell> spinor_basis_b {nuclear_framework_b, basisset_name};  // in non-orthogonal AO basis

    const auto K_a = spinor_basis_a.numberOfSpatialOrbitals();
    const auto K_b = spinor_basis_b.numberOfSpatialOrbitals();

    SquareMatrix<double> V_a = spinor_basis_a.quantize(Operator::NuclearAttraction(nuclear_framework_a)).parameters();
    SquareMatrix<double> V_b = spinor_basis_b.quantize(Operator::NuclearAttraction(nuclear_framework_b)).parameters();

    // T_a and T_b are equal to the corresponding block from the molecular kinetic integrals (T_a = T.block(0,0, K_a, K_a))
    SquareMatrix<double> T_a = spinor_basis_a.quantize(Operator::Kinetic()).parameters();
    SquareMatrix<double> T_b = spinor_basis_b.quantize(Operator::Kinetic()).parameters();

    // Create partition matrices for both atoms
    const auto p_a = SquareMatrix<double>::PartitionMatrix(0, K_a, K);
    const auto p_b = SquareMatrix<double>::PartitionMatrix(K_a, K_b, K);

    // Retrieve the molecular integrals
    const auto S = spinor_basis.quantize(Operator::Overlap()).parameters();
    const auto T = spinor_basis.quantize(Operator::Kinetic()).parameters();
    const auto V = spinor_basis.quantize(Operator::NuclearAttraction(molecule)).parameters();
    const auto g = spinor_basis.quantize(Operator::Coulomb()).parameters();

    SquareMatrix<double> H = T + V;

    // Decompose the integrals corresponding to the formula's in Mario's thesis
    SquareMatrix<double> h_a = SquareMatrix<double>::Zero(K);
    SquareMatrix<double> h_b = SquareMatrix<double>::Zero(K);

    h_a.block(0, 0, K_a, K_a) = T_a + V_a;
    h_b.block(K_a, K_a, K_b, K_b) = T_b + V_b;

    SquareMatrix<double> h_ab = H - h_a - h_b;

    auto g_a = g;
    auto g_b = g;
    auto g_ab = g;
    auto g_ba = g;

    g_a.contractWithMatrix<double>(p_a, 0);
    g_a.contractWithMatrix<double>(p_a, 2);

    g_b.contractWithMatrix<double>(p_b, 0);
    g_b.contractWithMatrix<double>(p_b, 2);

    g_ab.contractWithMatrix<double>(p_a, 0);
    g_ab.contractWithMatrix<double>(p_b, 2);

    g_ba.contractWithMatrix<double>(p_b, 0);
    g_ba.contractWithMatrix<double>(p_a, 2);

    QCRankFourTensor<double> g_abba = g_ab.Eigen() + g_ba.Eigen();

    RSQHamiltonian<double> HAA {ScalarRSQOneElectronOperator<double>(h_a), ScalarRSQTwoElectronOperator<double>(g_a)};
    RSQHamiltonian<double> HBB {ScalarRSQOneElectronOperator<double>(h_b), ScalarRSQTwoElectronOperator<double>(g_b)};
    RSQHamiltonian<double> HAB {ScalarRSQOneElectronOperator<double>(h_ab), ScalarRSQTwoElectronOperator<double>(g_abba)};
    RSQHamiltonian<double> HA {ScalarRSQOneElectronOperator<double>(h_a + h_ab / 2), ScalarRSQTwoElectronOperator<double>(g_a.Eigen() + 0.5 * g_abba.Eigen())};
    RSQHamiltonian<double> HB {ScalarRSQOneElectronOperator<double>(h_b + h_ab / 2), ScalarRSQTwoElectronOperator<double>(g_b.Eigen() + 0.5 * g_abba.Eigen())};

    std::vector<RSQHamiltonian<double>> net_atomic_parameters {HAA, HBB};
    std::vector<RSQHamiltonian<double>> interaction_parameters {HAB};
    std::vector<RSQHamiltonian<double>> atomic_parameters {HA, HB};

    return AtomicDecompositionParameters(RSQHamiltonian<double>::Molecular(spinor_basis, molecule),
                                         net_atomic_parameters, interaction_parameters, atomic_parameters);
}


}  // namespace GQCP
