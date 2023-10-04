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

#define BOOST_TEST_MODULE "NOCI_test"

#include <boost/test/unit_test.hpp>

#include "Basis/NonOrthogonalBasis/GNonOrthogonalStateBasis.hpp"
#include "Basis/SpinorBasis/GSpinorBasis.hpp"
#include "Basis/Transformations/GTransformation.hpp"
#include "Mathematical/Optimization/Eigenproblem/GeneralizedEigenproblemSolver.hpp"
#include "Molecule/Molecule.hpp"
#include "QCMethod/HF/GHF/GHF.hpp"
#include "QCMethod/HF/GHF/GHFSCFSolver.hpp"
#include "QCMethod/NOCI/NOCI.hpp"
#include "QCMethod/NOCI/NOCIEnvironment.hpp"
#include "QCModel/NOCI/NOCIExpansion.hpp"
#include "Utilities/complex.hpp"

/**
 *  Test the NOCI QC Method.
 */
BOOST_AUTO_TEST_CASE(NOCI) {
    // It was for H2, at 2.5au internuclear distance for the 6-31G basis set.
    const auto molecule = GQCP::Molecule::HChain(2, 2.5, 0);  // H2, 2.5 bohr apart.

    // The general spinor basis is also needed, as we require the overlap operator in AO basis.
    const GQCP::GSpinorBasis<double, GQCP::GTOShell> g_spinor_basis {molecule, "6-31G"};
    const auto S = g_spinor_basis.overlap();

    // Initialize some non-orthogonal "generalised states".
    GQCP::SquareMatrix<double> basis_state_1 {8};
    // clang-format off
    basis_state_1 << -0.07443693,  0.12036042, -0.13557067,  0.15517005,  0.13315100,  -0.03074946, -0.92997876, -0.93718779,
                     -0.07874922,  0.15086478, -0.68085546,  0.77423311,  0.08712373,   0.25266303,  0.848079  ,  0.89108911,
                     -0.24580188,  0.26338108,  0.09556297, -0.12178159,  0.91319299,   0.90475733, -0.03994767,  0.12839983,
                     -0.38944259,  0.4101685 ,  0.45214166, -0.58335985, -0.90125958,  -0.93270816, -0.16410814, -0.32074956,
                     -0.26338108, -0.24580188, -0.12178159, -0.09556297, -0.90475733,   0.91319299, -0.12839983, -0.03994767,
                     -0.4101685 , -0.38944259, -0.58335985, -0.45214166,  0.93270817,  -0.90125958,  0.32074956, -0.16410814,
                     -0.12036042, -0.07443693,  0.15517005,  0.13557067,  0.03074946,   0.13315101,  0.93718779, -0.92997876,
                     -0.15086478, -0.07874922,  0.77423311,  0.68085546, -0.25266303,   0.08712373, -0.89108911,  0.84807900;
    // clang-format on
    GQCP::SquareMatrix<double> basis_state_2 {8};
    // clang-format off
    basis_state_2 <<  0.25851329, -0.14539151, -0.17177142, -0.01126487,  0.81289875, -0.77260907,  0.50167389, -0.44422385,
                      0.36593356, -0.28669343, -0.84796858, -0.13503625, -0.62437698,  0.96771154, -0.55231929,  0.30317456,
                      0.25853403,  0.14539669,  0.17176599, -0.01126146,  0.81450567,  0.7709918 , -0.501289  , -0.44451308,
                      0.36597032,  0.28670189,  0.847938  , -0.13501526, -0.62639487, -0.96647128,  0.5520554 ,  0.30349133,
                      0.10076798, -0.23874662,  0.04823423,  0.17685836,  0.42013282, -0.48352714, -0.79642816,  0.8239984 ,
                      0.16561668, -0.35007843,  0.19502141,  0.90182842, -0.55545195,  0.39170258,  0.56753639, -0.94408827,
                     -0.10075937, -0.23872464,  0.0482368 , -0.17686313, -0.42104909, -0.4826058 , -0.79588057, -0.82460595,
                     -0.16560552, -0.35003836,  0.19503259, -0.9018579 ,  0.55619574,  0.39048771,  0.56690551,  0.94451894;
    // clang-format on
    GQCP::SquareMatrix<double> basis_state_3 {8};
    // clang-format off
    basis_state_3 <<  -0.265842  ,  0.17716735, -0.15969328, -0.00308706,  0.84741422, -0.81942223,  0.41366608, -0.36889812,
                      -0.36278694,  0.36406651, -0.80340861, -0.13144475, -0.65133121,  1.03324716, -0.44428872,  0.24605534,
                      -0.26584976, -0.17716927,  0.15969112, -0.00308558,  0.84933355,  0.81745234, -0.41355441, -0.36897416,
                      -0.36280035, -0.36406982,  0.80339638, -0.13143372, -0.65375477, -1.0317334 ,  0.4442138 ,  0.24613641,
                      -0.09736842,  0.22594293,  0.06676532,  0.17043252,  0.3439281 , -0.39253422, -0.84701679,  0.86052195,
                      -0.15038318,  0.32968947,  0.23916148,  0.91374992, -0.47007004,  0.32909412,  0.60131983, -0.98100354,
                       0.0973641 ,  0.22593416,  0.06676626, -0.17043417, -0.34482098, -0.39170837, -0.84682078, -0.86073623,
                       0.15037797,  0.32967347,  0.2391655 , -0.91376119,  0.47081908,  0.32796701,  0.60109455,  0.98115454;
    // clang-format on
    // Transform the matrices to the correct transformation type.
    const auto state1 = GQCP::GTransformation<double> {basis_state_1};
    const auto state2 = GQCP::GTransformation<double> {basis_state_2};
    const auto state3 = GQCP::GTransformation<double> {basis_state_3};

    // Create a vector out of these three basis states.
    std::vector<GQCP::GTransformation<double>> basis_vector {state1, state2, state3};

    // Create a non-orthogonal state basis, using the basis state vector, the overlap operator in AO basis and the number of occupied orbitals.
    const auto NOS_basis = GQCP::GNonOrthogonalStateBasis<double> {basis_vector, S, molecule.numberOfElectrons()};

    // Quantize the Hamiltonian in the general spinor basis.
    const auto sq_hamiltonian = g_spinor_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));

    // Create a dense solver and corresponding environment and put them together in the QCMethod.
    auto environment = GQCP::NOCIEnvironment::Dense(sq_hamiltonian, NOS_basis, molecule);
    auto solver = GQCP::GeneralizedEigenproblemSolver::Dense<double>();

    // We ask for three states (which is all of them in this example) to be found.
    const auto qc_structure = GQCP::QCMethod::NOCI<double, GQCP::GNonOrthogonalStateBasis<double>>(NOS_basis, 3).optimize(solver, environment);

    // Since we use three states as our basis, we will gain th ground state energy and two excited states.
    std::vector<double> reference_energies {-1.04749865, -0.92999724, -0.90565568};

    // Compare the reference energies with the calculated energies.
    BOOST_CHECK(std::abs(reference_energies[0] - qc_structure.groundStateEnergy()) < 1e-5);
    BOOST_CHECK(std::abs(reference_energies[1] - qc_structure.energy(1)) < 1e-5);
    BOOST_CHECK(std::abs(reference_energies[2] - qc_structure.energy(2)) < 1e-5);

    // We will also check the expansion coefficients for each of the calculated states.
    // First, we initialize a reference for each state.
    std::vector<double> reference_coefficients_1 {-0.55829124, 0.49701288, 0.00357357};
    std::vector<double> reference_coefficients_2 {1.5119597, 1.55598945, -0.02485513};
    std::vector<double> reference_coefficients_3 {-1.06431638, -11.22652471, 10.462225};

    // Compare the calculated coefficients to the reference.
    BOOST_CHECK(std::abs(reference_coefficients_1[0] - qc_structure.groundStateParameters().coefficient(0)) < 1e-5);
    BOOST_CHECK(std::abs(reference_coefficients_1[1] - qc_structure.groundStateParameters().coefficient(1)) < 1e-5);
    BOOST_CHECK(std::abs(reference_coefficients_1[2] - qc_structure.groundStateParameters().coefficient(2)) < 1e-5);
    BOOST_CHECK(std::abs(reference_coefficients_2[0] - qc_structure.parameters(1).coefficient(0)) < 1e-5);
    BOOST_CHECK(std::abs(reference_coefficients_2[1] - qc_structure.parameters(1).coefficient(1)) < 1e-5);
    BOOST_CHECK(std::abs(reference_coefficients_2[2] - qc_structure.parameters(1).coefficient(2)) < 1e-5);
    BOOST_CHECK(std::abs(reference_coefficients_3[0] - qc_structure.parameters(2).coefficient(0)) < 1e-5);
    BOOST_CHECK(std::abs(reference_coefficients_3[1] - qc_structure.parameters(2).coefficient(1)) < 1e-5);
    BOOST_CHECK(std::abs(reference_coefficients_3[2] - qc_structure.parameters(2).coefficient(2)) < 1e-5);
}


/**
 *  Test the NOCI model 1DM calculation.
 */
BOOST_AUTO_TEST_CASE(NOCI_1DM) {
    // It was for H2, at 2.5au internuclear distance for the 6-31G basis set.
    const auto molecule = GQCP::Molecule::HChain(2, 2.5, 0);  // H2, 2.5 bohr apart.

    // The general spinor basis is also needed, as we require the overlap operator in AO basis.
    const GQCP::GSpinorBasis<double, GQCP::GTOShell> g_spinor_basis {molecule, "6-31G"};
    const auto S = g_spinor_basis.overlap();

    // Initialize some non-orthogonal "generalised states".
    GQCP::SquareMatrix<double> basis_state_1 {8};
    // clang-format off
    basis_state_1 << -0.07443693,  0.12036042, -0.13557067,  0.15517005,  0.13315100,  -0.03074946, -0.92997876, -0.93718779,
                     -0.07874922,  0.15086478, -0.68085546,  0.77423311,  0.08712373,   0.25266303,  0.848079  ,  0.89108911,
                     -0.24580188,  0.26338108,  0.09556297, -0.12178159,  0.91319299,   0.90475733, -0.03994767,  0.12839983,
                     -0.38944259,  0.4101685 ,  0.45214166, -0.58335985, -0.90125958,  -0.93270816, -0.16410814, -0.32074956,
                     -0.26338108, -0.24580188, -0.12178159, -0.09556297, -0.90475733,   0.91319299, -0.12839983, -0.03994767,
                     -0.4101685 , -0.38944259, -0.58335985, -0.45214166,  0.93270817,  -0.90125958,  0.32074956, -0.16410814,
                     -0.12036042, -0.07443693,  0.15517005,  0.13557067,  0.03074946,   0.13315101,  0.93718779, -0.92997876,
                     -0.15086478, -0.07874922,  0.77423311,  0.68085546, -0.25266303,   0.08712373, -0.89108911,  0.84807900;
    // clang-format on
    GQCP::SquareMatrix<double> basis_state_2 {8};
    // clang-format off
    basis_state_2 <<  0.25851329, -0.14539151, -0.17177142, -0.01126487,  0.81289875, -0.77260907,  0.50167389, -0.44422385,
                      0.36593356, -0.28669343, -0.84796858, -0.13503625, -0.62437698,  0.96771154, -0.55231929,  0.30317456,
                      0.25853403,  0.14539669,  0.17176599, -0.01126146,  0.81450567,  0.7709918 , -0.501289  , -0.44451308,
                      0.36597032,  0.28670189,  0.847938  , -0.13501526, -0.62639487, -0.96647128,  0.5520554 ,  0.30349133,
                      0.10076798, -0.23874662,  0.04823423,  0.17685836,  0.42013282, -0.48352714, -0.79642816,  0.8239984 ,
                      0.16561668, -0.35007843,  0.19502141,  0.90182842, -0.55545195,  0.39170258,  0.56753639, -0.94408827,
                     -0.10075937, -0.23872464,  0.0482368 , -0.17686313, -0.42104909, -0.4826058 , -0.79588057, -0.82460595,
                     -0.16560552, -0.35003836,  0.19503259, -0.9018579 ,  0.55619574,  0.39048771,  0.56690551,  0.94451894;
    // clang-format on
    GQCP::SquareMatrix<double> basis_state_3 {8};
    // clang-format off
    basis_state_3 <<  -0.265842  ,  0.17716735, -0.15969328, -0.00308706,  0.84741422, -0.81942223,  0.41366608, -0.36889812,
                      -0.36278694,  0.36406651, -0.80340861, -0.13144475, -0.65133121,  1.03324716, -0.44428872,  0.24605534,
                      -0.26584976, -0.17716927,  0.15969112, -0.00308558,  0.84933355,  0.81745234, -0.41355441, -0.36897416,
                      -0.36280035, -0.36406982,  0.80339638, -0.13143372, -0.65375477, -1.0317334 ,  0.4442138 ,  0.24613641,
                      -0.09736842,  0.22594293,  0.06676532,  0.17043252,  0.3439281 , -0.39253422, -0.84701679,  0.86052195,
                      -0.15038318,  0.32968947,  0.23916148,  0.91374992, -0.47007004,  0.32909412,  0.60131983, -0.98100354,
                       0.0973641 ,  0.22593416,  0.06676626, -0.17043417, -0.34482098, -0.39170837, -0.84682078, -0.86073623,
                       0.15037797,  0.32967347,  0.2391655 , -0.91376119,  0.47081908,  0.32796701,  0.60109455,  0.98115454;
    // clang-format on
    // Transform the matrices to the correct transformation type.
    const auto state1 = GQCP::GTransformation<double> {basis_state_1};
    const auto state2 = GQCP::GTransformation<double> {basis_state_2};
    const auto state3 = GQCP::GTransformation<double> {basis_state_3};

    // Create a vector out of these three basis states.
    std::vector<GQCP::GTransformation<double>> basis_vector {state1, state2, state3};

    // Create a non-orthogonal state basis, using the basis state vector, the overlap operator in AO basis and the number of occupied orbitals.
    const auto NOS_basis = GQCP::GNonOrthogonalStateBasis<double> {basis_vector, S, molecule.numberOfElectrons()};

    // Quantize the Hamiltonian in the general spinor basis.
    const auto sq_hamiltonian = g_spinor_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));

    // Create a dense solver and corresponding environment and put them together in the QCMethod.
    auto environment = GQCP::NOCIEnvironment::Dense(sq_hamiltonian, NOS_basis, molecule);
    auto solver = GQCP::GeneralizedEigenproblemSolver::Dense<double>();

    // We do not specify the number of states, meaning we only request the ground state.
    const auto NOCI_parameters = GQCP::QCMethod::NOCI<double, GQCP::GNonOrthogonalStateBasis<double>>(NOS_basis).optimize(solver, environment).groundStateParameters();

    // Initialize a reference 1DM.
    GQCP::SquareMatrix<double> reference_1DM {8};
    // clang-format off
    reference_1DM <<  0.04286838,  0.05940042,  0.05857023,  0.08654403,  0.0287501 ,  0.04463992,  0.004545  ,  0.0040954 ,
                      0.05940042,  0.08403974,  0.07425935,  0.10822316,  0.05119187,  0.07919022,  0.01406753,  0.01633957,
                      0.05857023,  0.07425935,  0.10779418,  0.16531073, -0.00559086, -0.0075571 , -0.02407526, -0.035894  ,
                      0.08654403,  0.10822316,  0.16531073,  0.25449135, -0.01811907, -0.0262255 , -0.04223742, -0.06216742,
                      0.0287501 ,  0.05119187, -0.00559086, -0.01811907,  0.10258304,  0.15742264,  0.05787152,  0.07747057,
                      0.04463992,  0.07919022, -0.0075571 , -0.0262255 ,  0.15742264,  0.24160629,  0.08846938,  0.11835985,
                      0.004545  ,  0.01406753, -0.02407526, -0.04223742,  0.05787152,  0.08846938,  0.03797135,  0.05184621,
                      0.0040954 ,  0.01633957, -0.035894  , -0.06216742,  0.07747057,  0.11835985,  0.05184621,  0.0709642 ;
    // clang-format on
    // Check the calculated 1DM against the reference.
    BOOST_CHECK(reference_1DM.isApprox(NOCI_parameters.calculate1DM().matrix(), 1e-6));
}

/**
 *  Test the complex NOCI model 1DM calculation, restricted vs. unrestricted.
 */
BOOST_AUTO_TEST_CASE(NOCI_1DM_complex) {
    // It was for H2, at 2.5au internuclear distance for the 6-31G basis set.
    const auto molecule = GQCP::Molecule::HChain(2, 2.5, 0);  // H2, 2.5 bohr apart.

    // The restricted spin orbital basis is also needed, as we require the overlap operator in AO basis.
    const GQCP::RSpinOrbitalBasis<GQCP::complex, GQCP::GTOShell> r_spin_orbital_basis {molecule, "6-31G"};
    const auto r_S = r_spin_orbital_basis.overlap();

    // The unrestricted spin orbital basis is also needed, as we require the overlap operator in AO basis.
    const GQCP::USpinOrbitalBasis<GQCP::complex, GQCP::GTOShell> u_spin_orbital_basis {molecule, "6-31G"};
    const auto u_S = u_spin_orbital_basis.overlap();

    // Initialize two non-orthogonal "restricted states".
    GQCP::SquareMatrix<GQCP::complex> state_1 {4};
    using namespace std::complex_literals;
    // clang-format off
    state_1 << -0.07443693+0.0i,  0.12036042+0.0i, -0.13557067+0.0i,  0.15517005+0.0i,
               -0.07874922+0.0i,  0.15086478+0.0i, -0.68085546+0.0i,  0.77423311+0.0i,
               -0.24580188+0.0i,  0.26338108+0.0i,  0.09556297+0.0i, -0.12178159+0.0i,
               -0.38944259+0.0i,  0.4101685+0.0i ,  0.45214166+0.0i, -0.58335985+0.0i;
    // clang-format on
    GQCP::SquareMatrix<GQCP::complex> state_2 {4};
    // clang-format off
    state_2 <<  0.25851329+0.0i, -0.14539151+0.0i, -0.17177142+0.0i, -0.01126487+0.0i,
                0.36593356+0.0i, -0.28669343+0.0i, -0.84796858+0.0i, -0.13503625+0.0i,
                0.25853403+0.0i,  0.14539669+0.0i,  0.17176599+0.0i, -0.01126146+0.0i,
                0.36597032+0.0i,  0.28670189+0.0i,  0.847938+0.0i  , -0.13501526+0.0i;
    // clang-format on
    GQCP::SquareMatrix<GQCP::complex> state_3 {4};
    // clang-format off
    state_3 <<  -0.265842+0.0i  ,  0.17716735+0.0i, -0.15969328+0.0i, -0.00308706+0.0i,
                -0.36278694+0.0i,  0.36406651+0.0i, -0.80340861+0.0i, -0.13144475+0.0i,
                -0.26584976+0.0i, -0.17716927+0.0i,  0.15969112+0.0i, -0.00308558+0.0i,
                -0.36280035+0.0i, -0.36406982+0.0i,  0.80339638+0.0i, -0.13143372+0.0i;
    // clang-format on
    // Transform the matrices to the correct transformation type.
    const auto r_basis_state_1 = GQCP::RTransformation<GQCP::complex> {state_1};
    const auto r_basis_state_2 = GQCP::RTransformation<GQCP::complex> {state_2};
    const auto r_basis_state_3 = GQCP::RTransformation<GQCP::complex> {state_3};

    const auto u_basis_state_1 = GQCP::UTransformation<GQCP::complex>::FromRestricted(r_basis_state_1);
    const auto u_basis_state_2 = GQCP::UTransformation<GQCP::complex>::FromRestricted(r_basis_state_2);
    const auto u_basis_state_3 = GQCP::UTransformation<GQCP::complex>::FromRestricted(r_basis_state_3);

    // Create a vector out of these three basis states.
    std::vector<GQCP::RTransformation<GQCP::complex>> r_basis_vector {r_basis_state_1, r_basis_state_2, r_basis_state_3};
    std::vector<GQCP::UTransformation<GQCP::complex>> u_basis_vector {u_basis_state_1, u_basis_state_2, u_basis_state_3};

    // Create a non-orthogonal state basis, using the basis state vector, the overlap operator in AO basis and the number of occupied orbitals.
    const auto RNOS_basis = GQCP::RNonOrthogonalStateBasis<GQCP::complex> {r_basis_vector, r_S, molecule.numberOfElectronPairs()};
    const auto UNOS_basis = GQCP::UNonOrthogonalStateBasis<GQCP::complex> {u_basis_vector, u_S, molecule.numberOfElectronPairs(), molecule.numberOfElectronPairs()};

    // To evaluate the Hamiltonian in the non-orthogonal state basis, we need the second quantized Hamiltonian in AO basis.
    const auto rsq_hamiltonian = r_spin_orbital_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));
    const auto usq_hamiltonian = u_spin_orbital_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));

    // Create a dense solver and corresponding environment and put them together in the QCMethod.
    auto r_environment = GQCP::NOCIEnvironment::Dense(rsq_hamiltonian, RNOS_basis, molecule);
    auto u_environment = GQCP::NOCIEnvironment::Dense(usq_hamiltonian, UNOS_basis, molecule);
    auto solver = GQCP::GeneralizedEigenproblemSolver::Dense<GQCP::complex>();

    // We do not specify the number of states, meaning we only request the ground state.
    const auto RNOCI_parameters = GQCP::QCMethod::NOCI<GQCP::complex, GQCP::RNonOrthogonalStateBasis<GQCP::complex>>(RNOS_basis).optimize(solver, r_environment).groundStateParameters();
    const auto UNOCI_parameters = GQCP::QCMethod::NOCI<GQCP::complex, GQCP::UNonOrthogonalStateBasis<GQCP::complex>>(UNOS_basis).optimize(solver, u_environment).groundStateParameters();

    const auto orbitalRDM = RNOCI_parameters.calculate1DM().matrix();
    const auto resolvedRDM_a = UNOCI_parameters.calculate1DM().alpha().matrix();
    const auto resolvedRDM_b = UNOCI_parameters.calculate1DM().beta().matrix();

    // The restricted 1DM should equal the alpha/beta 1DM of the unrestricted case.
    BOOST_CHECK(orbitalRDM.isApprox(resolvedRDM_a, 1e-6));
    BOOST_CHECK(orbitalRDM.isApprox(resolvedRDM_b, 1e-6));
}


/**
 *  Test the NOCI QC Method.
 */
BOOST_AUTO_TEST_CASE(NOCI_complex_generalized) {
    // This is a self implemented test case, with rotated H3 states in STO-3G.
    const auto molecule = GQCP::Molecule::HRingFromDistance(3, 1.88973, 0);  // H3, 1 Angstrom apart.

    // The general spinor basis is also needed, as we require the overlap operator in AO basis.
    const GQCP::GSpinorBasis<GQCP::complex, GQCP::GTOShell> g_spinor_basis {molecule, "STO-3G"};
    const auto S = g_spinor_basis.overlap();

    // Initialize some non-orthogonal "generalised states".
    GQCP::SquareMatrix<GQCP::complex> basis_state_1 {6};
    using namespace std::complex_literals;
    // clang-format off
    basis_state_1 << 0.460055771  - 5.63405828e-17i,  1.35013939  + 0.0i           ,  0.996501939    - 1.22036291e-16i,  0.535358621 - 6.55625222e-17i,  1.35013939 + 0.0i            ,  0.0 + 0.0i,
                     0.460055767  - 5.63405822e-17i,  1.35013943  + 0.0i           , -0.996501918    + 1.22036288e-16i,  0.535358665 - 6.55625275e-17i,  1.35013943 + 0.0i            ,  0.0 + 0.0i,
                     0.301611955  - 3.69368116e-17i,  3.07322180  + 0.0i           , -2.63772351e-08 + 3.23027965e-24i, -1.18334547  + 1.44918025e-16i,  3.07322180 + 0.0i            ,  0.0 + 0.0i,
                    -0.0956361988 + 0.0i           , -0.280666403 - 3.43717213e-17i, -0.207152401    + 0.0i           , -0.111290123 + 0.0i           , -0.280666403 - 3.43717213e-17i,  0.0 + 0.0i,
                    -0.0956361979 + 0.0i           , -0.280666413 - 3.43717224e-17i,  0.207152397    + 0.0i           , -0.111290132 + 0.0i           , -0.280666413 - 3.43717224e-17i,  0.0 + 0.0i,
                    -0.0626989655 + 0.0i           , -0.638860046 - 7.82377910e-17i,  5.48328845e-09 + 0.0i           ,  0.245993356 + 0.0i           , -0.638860046 - 7.82377910e-17i,  0.0 + 0.0i;
    // clang-format on
    GQCP::SquareMatrix<GQCP::complex> basis_state_2 {6};
    // clang-format off
    basis_state_2 <<  0.230027886 - 0.398419985i,  0.799802113 + 0.0i        ,  0.498250970    - 0.862995994i   ,  0.267679311 - 0.463634166i,  0.799802113 + 0.0i        , 0.0 + 0.0i,
                      0.230027883 - 0.398419981i,  0.799802140 + 0.0i        , -0.498250959    + 0.862995976i   ,  0.267679332 - 0.463634204i,  0.799802140 + 0.0i        , 0.0 + 0.0i,
                      0.150805978 - 0.261203615i,  1.82053003  + 0.0i        , -1.31886175e-08 + 2.28433557e-08i, -0.591672737 + 1.02480724i ,  1.82053003  + 0.0i        , 0.0 + 0.0i,
                     -0.161442683 + 0.0i        , -0.140333202 - 0.243064235i, -0.349692268    + 0.0i           , -0.187867944 + 0.0i        , -0.140333202 - 0.243064235i, 0.0 + 0.0i,
                     -0.161442681 + 0.0i        , -0.140333206 - 0.243064243i,  0.349692261    + 0.0i           , -0.187867959 + 0.0i        , -0.140333206 - 0.243064243i, 0.0 + 0.0i,
                     -0.105841609 + 0.0i        , -0.319430023 - 0.553269029i,  9.25629424e-09 + 0.0i           ,  0.415259365 + 0.0i        , -0.319430023 - 0.553269029i, 0.0 + 0.0i;
    // clang-format on
    GQCP::SquareMatrix<GQCP::complex> basis_state_3 {6};
    // clang-format off
    basis_state_3 <<  -1.02152902e-16 - 0.460055771i,  0.615580024    + 0.0i        , -2.21267879e-16 - 0.996501939i   , -1.18873494e-16 - 0.535358621i, 0.615580024    +0.0i         ,  0.0 + 0.0i,
                      -1.02152901e-16 - 0.460055767i,  0.615580044    + 0.0i        ,  2.21267875e-16 + 0.996501918i   , -1.18873503e-16 - 0.535358665i, 0.615580044    +0.0i         ,  0.0 + 0.0i,
                      -6.69713075e-17 - 0.301611955i,  1.40119899     + 0.0i        ,  5.85692274e-24 + 2.63772351e-08i,  2.62755478e-16 + 1.18334547i , 1.40119899     +0.0i         ,  0.0 + 0.0i,
                      -0.209756967    + 0.0i        ,  6.23204607e-17 - 0.280666403i, -0.454343228    + 0.0i           , -0.244090407    +0.0i         , 6.23204607e-17 - 0.280666403i,  0.0 + 0.0i,
                      -0.209756965    + 0.0i        ,  6.23204627e-17 - 0.280666413i,  0.454343219    + 0.0i           , -0.244090427    +0.0i         , 6.23204627e-17 - 0.280666413i,  0.0 + 0.0i,
                      -0.137516390    + 0.0i        ,  1.41855426e-16 - 0.638860046i,  1.20263872e-08 + 0.0i           ,  0.539532320    +0.0i         , 1.41855426e-16 - 0.638860046i,  0.0 + 0.0i;
    // clang-format on
    // Transform the matrices to the correct transformation type.
    const auto state1 = GQCP::GTransformation<GQCP::complex> {basis_state_1};
    const auto state2 = GQCP::GTransformation<GQCP::complex> {basis_state_2};
    const auto state3 = GQCP::GTransformation<GQCP::complex> {basis_state_3};

    // Create a vector out of these three basis states.
    std::vector<GQCP::GTransformation<GQCP::complex>> basis_vector {state1, state2, state3};

    // Create a non-orthogonal state basis, using the basis state vector, the overlap operator in AO basis and the number of occupied orbitals.
    const auto NOS_basis = GQCP::GNonOrthogonalStateBasis<GQCP::complex> {basis_vector, S, molecule.numberOfElectrons()};

    // Quantize the Hamiltonian in the general spinor basis.
    const auto sq_hamiltonian = g_spinor_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));

    // Create a dense solver and corresponding environment and put them together in the QCMethod.
    auto environment = GQCP::NOCIEnvironment::Dense(sq_hamiltonian, NOS_basis, molecule);
    auto solver = GQCP::GeneralizedEigenproblemSolver::Dense<GQCP::complex>();

    // We ask for three states (which is all of them in this example) to be found.
    const auto qc_structure = GQCP::QCMethod::NOCI<GQCP::complex, GQCP::GNonOrthogonalStateBasis<GQCP::complex>>(NOS_basis).optimize(solver, environment);

    // Initialize the reference energy from a non-related diagonalization of the constructed Hamiltonian.
    const double reference = -0.97341735;

    // Check the result with the reference.
    BOOST_CHECK(std::abs(reference - qc_structure.groundStateEnergy()) < 1e-6);
}

/**
 *  Test the UNOCI in the case where two zero overlap values are present.
 */
BOOST_AUTO_TEST_CASE(NOCI_unrestricted_two_zero) {
    // It was for H3, at 1Å internuclear distance for the STO-3G basis set.
    const auto molecule = GQCP::Molecule::HRingFromDistance(3, 1.8897259886);  // H3, 1Å apart.

    // The unrestricted spin orbital basis is also needed, as we require the overlap operator in AO basis.
    const GQCP::USpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "STO-3G"};
    const auto S = spin_orbital_basis.overlap();

    // Initialize two non-orthogonal "unrestricted states".
    GQCP::SquareMatrix<double> state_1_a {3};
    GQCP::SquareMatrix<double> state_1_b {3};
    // clang-format off
    state_1_a << 0.301612395,  0.0        ,  1.18334666 ,
                 0.460055205, -0.996503166, -0.535359699,
                 0.460055205,  0.996503166, -0.535359699;

    state_1_b << 0.638858962, -1.04073944 ,  0.0        ,
                 0.280666711,  0.647678158, -0.996503166,
                 0.280666711,  0.647678158,  0.996503166;
    // clang-format on
    GQCP::SquareMatrix<double> state_2_a {3};
    GQCP::SquareMatrix<double> state_2_b {3};
    // clang-format off
    state_2_a <<  0.30161239,  0.63885896,  1.18334666,
                  0.46005521,  0.28066671, -0.5353597 ,
                  0.46005521,  0.28066671, -0.5353597 ;

    state_2_b <<  0.0        , -1.04073944 ,  0.0        ,
                 -0.996503166,  0.647678158, -0.996503166,
                  0.996503166,  0.647678158,  0.996503166;
    // clang-format on

    // Transform the matrices to the correct transformation type.
    const auto basis_state_1 = GQCP::UTransformation<double> {GQCP::UTransformationComponent<double> {state_1_a}, GQCP::UTransformationComponent<double> {state_1_b}};
    const auto basis_state_2 = GQCP::UTransformation<double> {GQCP::UTransformationComponent<double> {state_2_a}, GQCP::UTransformationComponent<double> {state_2_b}};

    // Create a vector out of these three basis states.
    std::vector<GQCP::UTransformation<double>> basis_vector {basis_state_1, basis_state_2};

    // Create a non-orthogonal state basis, using the basis state vector, the overlap operator in AO basis and the number of occupied orbitals for alpha and beta.
    const auto NOS_basis = GQCP::UNonOrthogonalStateBasis<double> {basis_vector, S, 2, 1};

    // To evaluate the Hamiltonian in the non-orthogonal state basis, we need the second quantized Hamiltonian in AO basis.
    const auto sq_hamiltonian = spin_orbital_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));

    // Create a dense solver and corresponding environment and put them together in the QCMethod.
    auto environment = GQCP::NOCIEnvironment::Dense(sq_hamiltonian, NOS_basis, molecule);
    auto solver = GQCP::GeneralizedEigenproblemSolver::Dense<double>();

    // We do not specify the number of states, meaning we only request the ground state.
    const auto NOCI_model = GQCP::QCMethod::NOCI<double, GQCP::UNonOrthogonalStateBasis<double>>(NOS_basis).optimize(solver, environment);

    // Initialize the reference energy from a non-related diagonalization of the constructed Hamiltonian.
    const auto reference_energy = -1.34344577;

    // Check the energy versus the reference.
    BOOST_CHECK(std::abs(reference_energy - NOCI_model.groundStateEnergy()) < 1e-6);
}

/**
 *  Test the UNOCI in the case where one zero overlap value is present.
 */
BOOST_AUTO_TEST_CASE(NOCI_unrestricted_one_zero) {
    // It was for H3, at 1.6Å internuclear distance for the STO-3G basis set.
    const auto molecule = GQCP::Molecule::HRingFromDistance(3, 3.023561581760001);  // H3, 1.6Å apart.

    // The unrestricted spin orbital basis is also needed, as we require the overlap operator in AO basis.
    const GQCP::USpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "STO-3G"};
    const auto S = spin_orbital_basis.overlap();

    // Initialize two non-orthogonal "unrestricted states".
    GQCP::SquareMatrix<double> state_1 {3};
    // clang-format off
    state_1 << 0.639638709, -0.801711730,  0.189534772,
               0.0        ,  0.0        , -1.04297785 ,
               0.639638715,  0.801711722,  0.189534786;
    // clang-format on
    GQCP::SquareMatrix<double> state_2 {3};
    // clang-format off
    state_2 <<  0.0,  0.667128998, -0.801711735,
                1.0, -0.296315367,  0.0        ,
                0.0,  0.667129020,  0.801711717;
    // clang-format on

    // Transform the matrices to the correct transformation type.
    const auto basis_state_1 = GQCP::UTransformation<double> {GQCP::UTransformationComponent<double> {state_1}, GQCP::UTransformationComponent<double> {state_2}};
    const auto basis_state_2 = GQCP::UTransformation<double> {GQCP::UTransformationComponent<double> {state_2}, GQCP::UTransformationComponent<double> {state_1}};

    // Create a vector out of these three basis states.
    std::vector<GQCP::UTransformation<double>> basis_vector {basis_state_1, basis_state_2};

    // Create a non-orthogonal state basis, using the basis state vector, the overlap operator in AO basis and the number of occupied orbitals for alpha and beta.
    const auto NOS_basis = GQCP::UNonOrthogonalStateBasis<double> {basis_vector, S, 2, 1};

    // To evaluate the Hamiltonian in the non-orthogonal state basis, we need the second quantized Hamiltonian in AO basis.
    const auto sq_hamiltonian = spin_orbital_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));

    // Create a dense solver and corresponding environment and put them together in the QCMethod.
    auto environment = GQCP::NOCIEnvironment::Dense(sq_hamiltonian, NOS_basis, molecule);
    auto solver = GQCP::GeneralizedEigenproblemSolver::Dense<double>();

    // We do not specify the number of states, meaning we only request the ground state.
    const auto NOCI_model = GQCP::QCMethod::NOCI<double, GQCP::UNonOrthogonalStateBasis<double>>(NOS_basis).optimize(solver, environment);

    // Initialize the reference energy from a non-related diagonalization of the constructed Hamiltonian.
    const auto reference_energy = -1.37371;

    // Check the energy versus the reference.
    BOOST_CHECK(std::abs(reference_energy - NOCI_model.groundStateEnergy()) < 1e-6);
}

/**
 *  This test checks whether the lower lying complex GHF solution can indeed be found.
 *  Note that this solution can also be found using real valued parameters.
 */
BOOST_AUTO_TEST_CASE(h3_sto3g_s2) {

    // Do our own GHF calculation.
    const auto molecule = GQCP::Molecule::HRingFromDistance(3, 1.8897259886);  // H3-ring, 1 Angstrom apart.

    const GQCP::GSpinorBasis<double, GQCP::GTOShell> spinor_basis {molecule, "sto-3g"};
    const auto sq_hamiltonian = spinor_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));  // In an AO basis.
    const auto S = spinor_basis.overlap();

    // Solve the GHF SCF equations, using a special initial guess from @xdvriend's implementation. This makes sure the off-diagonal spin-blocks are not zero blocks and helps the calculation to converge to a true GHF solution.
    GQCP::SquareMatrix<double> C_initial_matrix {6};
    // clang-format off
    C_initial_matrix << -0.3100721,  -0.15761163, -0.51612194, -0.38100148,  0.57090929, -0.37620802,
                        -0.00741269,  0.38801568, -0.25974834, -0.41043789, -0.67141074, -0.40332126,
                        -0.61961507,  0.18043708,  0.58367365,  0.17317687,  0.05464039, -0.45811451,
                         0.67031756,  0.28266352,  0.37079814, -0.23639173,  0.37758712, -0.3671939,
                         0.18059725, -0.8326703,   0.16282789, -0.03436191, -0.27832567, -0.41095738,
                         0.19477298,  0.13713633, -0.4018331,   0.77416187,  0.01572939, -0.42686445;
    // clang-format on

    const GQCP::GTransformation<double> C_initial {C_initial_matrix};
    GQCP::GHFSCFEnvironment<double> environment {3, sq_hamiltonian, S, C_initial};

    auto solver = GQCP::GHFSCFSolver<double>::Plain(1.0e-08, 4000);
    const auto qc_structure = GQCP::QCMethod::GHF<double>().optimize(solver, environment);
    const auto nuc_rep = GQCP::NuclearRepulsionOperator(molecule.nuclearFramework()).value();

    const auto S2 = spinor_basis.quantize(GQCP::ElectronicSpinSquaredOperator());
    const auto Sz = spinor_basis.quantize(GQCP::ElectronicSpin_zOperator());

    const auto D_ghf = qc_structure.groundStateParameters().calculateScalarBasis1DM();
    const auto d_ghf = qc_structure.groundStateParameters().calculateScalarBasis2DM();

    // Create a non-orthogonal state basis, using the basis state vector, the overlap operator in AO basis and the number of occupied orbitals for alpha and beta.
    std::vector<GQCP::GTransformation<double>> basis_vector {qc_structure.groundStateParameters().expansion()};
    const auto NOS_basis = GQCP::GNonOrthogonalStateBasis<double> {basis_vector, S, 3};

    // Create a dense solver and corresponding environment and put them together in the QCMethod.
    auto noci_environment = GQCP::NOCIEnvironment::Dense(sq_hamiltonian, NOS_basis, molecule);
    auto noci_solver = GQCP::GeneralizedEigenproblemSolver::Dense<double>();

    // We do not specify the number of states, meaning we only request the ground state.
    const auto NOCI_model = GQCP::QCMethod::NOCI<double, GQCP::GNonOrthogonalStateBasis<double>>(NOS_basis).optimize(noci_solver, noci_environment);

    const auto D_noci = NOCI_model.groundStateParameters().calculate1DM();
    const auto d_noci = NOCI_model.groundStateParameters().calculate2DM();

    // Initialize a reference energy. (From the code of @xdvriend.)
    const double reference_energy = -1.34044;

    // Check if the converged energy matches the reference energy.
    BOOST_CHECK(std::abs(reference_energy - NOCI_model.groundStateEnergy()) < 1e-6);
    BOOST_CHECK(d_ghf.tensor().isApprox(d_noci.tensor()));
    BOOST_CHECK(std::abs(0.840668 - S2.calculateExpectationValue(D_noci, d_noci)) < 1e-6);
}
