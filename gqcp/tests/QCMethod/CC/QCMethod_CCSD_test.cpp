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

#define BOOST_TEST_MODULE "CCSD"

#include <boost/test/unit_test.hpp>

#include "Basis/Transformations/transform.hpp"
#include "ONVBasis/SpinUnresolvedONV.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/CC/CCSD.hpp"
#include "QCMethod/CC/CCSDEnvironment.hpp"
#include "QCMethod/CC/CCSDSolver.hpp"
#include "QCMethod/HF/RHF/DiagonalRHFFockMatrixObjective.hpp"
#include "QCMethod/HF/RHF/RHF.hpp"
#include "QCMethod/HF/RHF/RHFSCFSolver.hpp"


/**
 *  Check if the implementation of spinor-CCSD is correct, by comparing with a reference by crawdad (https://github.com/CrawfordGroup/ProgrammingProjects/tree/master/Project%2305).
 *
 *  The system under consideration is H2O in an STO-3G basisset.
 */
BOOST_AUTO_TEST_CASE(h2o_crawdad) {

    // Prepare the canonical RHF spin-orbital basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o_crawdad.xyz");
    const auto N = molecule.numberOfElectrons();

    GQCP::RSpinorBasis<double, GQCP::GTOShell> r_spinor_basis {molecule, "STO-3G"};
    const auto r_sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(r_spinor_basis, molecule);  // in an AO basis
    const auto K = r_spinor_basis.numberOfSpatialOrbitals();

    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(N, r_sq_hamiltonian, r_spinor_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective {r_sq_hamiltonian};
    const auto rhf_qc_structure = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment);
    const auto rhf_parameters = rhf_qc_structure.groundStateParameters();

    r_spinor_basis.transform(rhf_parameters.coefficientMatrix());


    // Check if the intermediate RHF results are correct. We can't continue if this isn't the case.
    const auto rhf_energy = rhf_qc_structure.groundStateEnergy() + GQCP::Operator::NuclearRepulsion(molecule).value();
    const double ref_rhf_energy = -74.942079928192;
    BOOST_REQUIRE(std::abs(rhf_energy - ref_rhf_energy) < 1.0e-09);


    // Create a GSpinorBasis since we have implement spinor-CCSD, and quantize the molecular Hamiltonian in it.
    const auto g_spinor_basis = GQCP::GSpinorBasis<double, GQCP::GTOShell>::FromRestricted(r_spinor_basis);
    const auto M = g_spinor_basis.numberOfSpinors();
    const auto g_sq_hamiltonian = GQCP::GSQHamiltonian<double>::Molecular(g_spinor_basis, molecule);  // in the canonical restricted spin-orbitals

    // Create the GHF ONV (which is actually just the RHF ONV, since we're using the canonical RHF orbitals) and the corresponding orbital space.
    const auto reference_onv = GQCP::SpinUnresolvedONV::GHF(2 * K, N, rhf_parameters.spinOrbitalEnergiesBlocked());
    const auto orbital_space = reference_onv.orbitalSpace();


    // Initialize an environment suitable for CCSD.
    auto environment = GQCP::CCSDEnvironment<double>::PerturbativeCCSD(g_sq_hamiltonian, orbital_space);

    // Since we're working with a Hartree-Fock reference, the perturbative amplitudes actually correspond to the MP2 amplitudes. This means that the initial CCSD energy correction is the MP2 energy correction.
    const double ref_mp2_correction_energy = -0.049149636120;
    const auto& t1 = environment.t1_amplitudes.back();

    BOOST_REQUIRE(t1.asImplicitMatrixSlice().asMatrix().isZero(1.0e-08));  // for a HF reference, the perturbative T1 amplitudes are zero

    const auto initial_ccsd_correction_energy = environment.electronic_energies.back();
    BOOST_REQUIRE(std::abs(initial_ccsd_correction_energy - ref_mp2_correction_energy) < 1.0e-10);


    // Prepare the CCSD solver and optimize the CCSD model parameters.
    auto solver = GQCP::CCSDSolver<double>::Plain();
    const auto ccsd_qc_structure = GQCP::QCMethod::CCSD<double>().optimize(solver, environment);

    const auto ccsd_correlation_energy = ccsd_qc_structure.groundStateEnergy();

    const double ref_ccsd_correlation_energy = -0.070680088376;
    BOOST_CHECK(std::abs(ccsd_correlation_energy - ref_ccsd_correlation_energy) < 1.0e-08);
}
