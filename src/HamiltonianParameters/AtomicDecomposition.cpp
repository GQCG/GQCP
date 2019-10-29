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

#include "Basis/SpinorBasis/RSpinorBasis.hpp"
#include "Operator/FirstQuantized/Operator.hpp"


namespace GQCP {


/**
 *  @param molecular_hamiltonian_parameters     the complete molecular Hamiltonian
 *  @param net_atomic_parameters                collection of the net atomic Hamiltonian
 *  @param interaction_parameters               collection of the atomic interaction Hamiltonian
 *  @param atomic_parameters                    collection of the atomic Hamiltonian
 */
AtomicDecompositionParameters::AtomicDecompositionParameters (const SQHamiltonian<double>& molecular_hamiltonian_parameters, const std::vector<SQHamiltonian<double>>& net_atomic_parameters, const std::vector<SQHamiltonian<double>>& interaction_parameters, const std::vector<SQHamiltonian<double>>& atomic_parameters) :
    molecular_hamiltonian_parameters (molecular_hamiltonian_parameters),
    net_atomic_parameters (net_atomic_parameters),
    interaction_parameters (interaction_parameters),
    atomic_parameters (atomic_parameters)
{}


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

    const RSpinorBasis<double, GTOShell> sp_basis (molecule, basisset_name);
    const auto K = sp_basis.numberOfSpatialOrbitals();


    // Retrieve an AO basis for the individual atoms so that we can retrieve net atomic nuclear integrals
    const NuclearFramework nuclear_framework_a ({atoms[0]});
    const NuclearFramework nuclear_framework_b ({atoms[1]});

    RSpinorBasis<double, GTOShell> sp_basis_a (nuclear_framework_a, basisset_name);  // in non-orthogonal AO basis
    RSpinorBasis<double, GTOShell> sp_basis_b (nuclear_framework_b, basisset_name);  // in non-orthogonal AO basis

    const auto K_a = sp_basis_a.numberOfSpatialOrbitals();
    const auto K_b = sp_basis_b.numberOfSpatialOrbitals();

    QCMatrix<double> V_a = sp_basis_a.quantize(Operator::NuclearAttraction(nuclear_framework_a)).parameters();
    QCMatrix<double> V_b = sp_basis_b.quantize(Operator::NuclearAttraction(nuclear_framework_b)).parameters();

    // T_a and T_b are equal to the corresponding block from the molecular kinetic integrals (T_a = T.block(0,0, K_a, K_a))
    QCMatrix<double> T_a = sp_basis_a.quantize(Operator::Kinetic()).parameters();
    QCMatrix<double> T_b = sp_basis_b.quantize(Operator::Kinetic()).parameters();

    // Create partition matrices for both atoms
    const auto p_a = SquareMatrix<double>::PartitionMatrix(0, K_a, K);
    const auto p_b = SquareMatrix<double>::PartitionMatrix(K_a, K_b, K);

    // Retrieve the molecular integrals
    const auto S = sp_basis.quantize(Operator::Overlap()).parameters();
    const auto T = sp_basis.quantize(Operator::Kinetic()).parameters();
    const auto V = sp_basis.quantize(Operator::NuclearAttraction(molecule)).parameters();
    const auto g = sp_basis.quantize(Operator::Coulomb()).parameters();

    QCMatrix<double> H = T + V;

    // Decompose the integrals corresponding to the formula's in Mario's thesis
    QCMatrix<double> h_a = QCMatrix<double>::Zero(K, K);
    QCMatrix<double> h_b = QCMatrix<double>::Zero(K, K);

    h_a.block(0, 0, K_a, K_a) = T_a + V_a;
    h_b.block(K_a , K_a, K_b, K_b) = T_b + V_b;

    QCMatrix<double> h_ab = H - h_a - h_b;

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

    QCRankFourTensor<double> g_abba = g_ab.Eigen() + g_ba.Eigen();

    SQHamiltonian<double> HAA (ScalarSQOneElectronOperator<double>({h_a}), ScalarSQTwoElectronOperator<double>({g_a}));
    SQHamiltonian<double> HBB (ScalarSQOneElectronOperator<double>({h_b}), ScalarSQTwoElectronOperator<double>({g_b}));
    SQHamiltonian<double> HAB (ScalarSQOneElectronOperator<double>({h_ab}), ScalarSQTwoElectronOperator<double>({g_abba}));
    SQHamiltonian<double> HA (ScalarSQOneElectronOperator<double>({h_a + h_ab/2}), ScalarSQTwoElectronOperator<double>({g_a.Eigen() + (0.5)*g_abba.Eigen()}));
    SQHamiltonian<double> HB (ScalarSQOneElectronOperator<double>({h_b + h_ab/2}), ScalarSQTwoElectronOperator<double>({g_b.Eigen() + (0.5)*g_abba.Eigen()}));

    std::vector<SQHamiltonian<double>> net_atomic_parameters = {HAA, HBB};
    std::vector<SQHamiltonian<double>> interaction_parameters = {HAB};
    std::vector<SQHamiltonian<double>> atomic_parameters = {HA, HB};

    return AtomicDecompositionParameters(SQHamiltonian<double>::Molecular(sp_basis, molecule), net_atomic_parameters, interaction_parameters, atomic_parameters);
}


}  // namespace GQCP
