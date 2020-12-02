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

#define BOOST_TEST_MODULE "SpinResolved1DM"

#include <boost/test/unit_test.hpp>

#include "DensityMatrix/SpinResolved1DM.hpp"
#include "Mathematical/Optimization/Eigenproblem/EigenproblemSolver.hpp"
#include "ONVBasis/SpinResolvedONVBasis.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/CI/CI.hpp"
#include "QCMethod/CI/CIEnvironment.hpp"
#include "QCMethod/HF/UHF/UHF.hpp"
#include "QCMethod/HF/UHF/UHFSCFSolver.hpp"
#include "QCModel/HF/UHF.hpp"


/**
 *  Check if the trace of the spin-resolved 1-DM yields the appropriate number of electrons.
 *
 *  The system of interested is H2O//STO-3G, with 7 spatial orbitals and a Fock space dimension of 441.
 */
BOOST_AUTO_TEST_CASE(trace) {

    // Set up the molecular Hamiltonian in a Löwdin-orthonormalized spinor basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o_Psi4_GAMESS.xyz");
    const size_t N_alpha = molecule.numberOfElectronPairs();
    const size_t N_beta = molecule.numberOfElectronPairs();
    const auto N = molecule.numberOfElectrons();

    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    const auto K = spinor_basis.numberOfSpatialOrbitals();
    spinor_basis.lowdinOrthonormalize();

    const auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spinor_basis, molecule);


    // Do a dense FCI calculation.
    const GQCP::SpinResolvedONVBasis onv_basis {K, N_alpha, N_beta};

    auto environment = GQCP::CIEnvironment::Dense(sq_hamiltonian, onv_basis);
    auto solver = GQCP::EigenproblemSolver::Dense();

    const auto linear_expansion = GQCP::QCMethod::CI<GQCP::SpinResolvedONVBasis>(onv_basis).optimize(solver, environment).groundStateParameters();


    // Calculate the 1-DMs, calculate the traces and check if they match the expected result.
    const auto D = linear_expansion.calculateSpinResolved1DM();

    BOOST_CHECK(std::abs(D.alpha().trace() - N_alpha) < 1.0e-12);
    BOOST_CHECK(std::abs(D.beta().trace() - N_beta) < 1.0e-12);
    BOOST_CHECK(std::abs(D.orbitalDensity().trace() - (N_alpha + N_beta)) < 1.0e-12);
}


/**
 *  Test if the expectation value of a one-electron operator in different orbital bases is the same.
 */
BOOST_AUTO_TEST_CASE(one_electron_operator_expectation_value_different_orbital_bases) {

    // Prepare the molecular Hamiltonian in the AO basis, in order to proceed with an UHF SCF calculation.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/ch4_crawdad.xyz");
    const auto N_a = molecule.numberOfElectronPairs();
    const auto N_b = molecule.numberOfElectronPairs();

    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "STO-3G"};
    const auto S = spin_orbital_basis.overlap();

    const auto hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spin_orbital_basis, molecule);  // In the AO basis.
    const auto K = hamiltonian.numberOfOrbitals();

    // Do the UHF SCF calculation to retrieve the UHF MOs.
    auto uhf_environment = GQCP::UHFSCFEnvironment<double>::WithCoreGuess(N_a, N_b, hamiltonian, S.parameters());
    auto diis_uhf_scf_solver = GQCP::UHFSCFSolver<double>::DIIS();

    const auto uhf_parameters = GQCP::QCMethod::UHF<double>().optimize(diis_uhf_scf_solver, uhf_environment).groundStateParameters();


    // Prepare three two-electron operators in different orbital bases.
    auto u_spin_orbital_basis = GQCP::USpinOrbitalBasis<double, GQCP::GTOShell>::FromRestricted(spin_orbital_basis);
    u_spin_orbital_basis.transform(uhf_parameters.expansion());
    auto u_hamiltonian = GQCP::USQHamiltonian<double>::Molecular(u_spin_orbital_basis, molecule);

    const auto& h_core_AO = u_hamiltonian.core();
    const auto h_core_MO = h_core_AO.transformed(uhf_parameters.expansion());

    const auto T_random = GQCP::UTransformation<double>::Random(K);
    const auto h_core_random = h_core_AO.transformed(T_random);

    // Prepare three density matrices in the corresponding orbital bases.
    const auto D_AO = uhf_parameters.calculateOrthonormalBasis1DM();
    const auto D_MO = D_AO.transformed(uhf_parameters.expansion());
    const auto D_random = D_AO.transformed(T_random);


    // Check if the expectation values match.
    const double exp_val_AO = h_core_AO.calculateExpectationValue(D_AO);
    const double exp_val_MO = h_core_MO.calculateExpectationValue(D_MO);
    const double exp_val_random = h_core_random.calculateExpectationValue(D_random);

    BOOST_CHECK(std::abs(exp_val_AO - exp_val_MO) < 1.0e-12);
    BOOST_CHECK(std::abs(exp_val_AO - exp_val_random) < 1.0e-12);
}
