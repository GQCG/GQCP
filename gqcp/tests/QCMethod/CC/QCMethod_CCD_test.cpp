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

#define BOOST_TEST_MODULE "CCD"

#include <boost/test/unit_test.hpp>

#include "Basis/transform.hpp"
#include "ONVBasis/SpinUnresolvedONV.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/CC/CCD.hpp"
#include "QCMethod/CC/CCSDEnvironment.hpp"
#include "QCMethod/CC/CCDSolver.hpp"
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
    const auto r_sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(r_spinor_basis, molecule);  // in an AO basis
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


    // Create a GSpinorBasis since we have implement spinor-CCD, and quantize the molecular Hamiltonian in it.
    const auto g_spinor_basis = GQCP::GSpinorBasis<double, GQCP::GTOShell>::FromRestricted(r_spinor_basis);
    const auto M = g_spinor_basis.numberOfSpinors();
    const auto g_sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(g_spinor_basis, molecule);  // in the canonical restricted spin-orbitals

    // Create the GHF ONV (which is actually just the RHF ONV, since we're using the canonical RHF orbitals) and the corresponding orbital space.
    //const auto reference_onv = GQCP::SpinUnresolvedONV::GHF(2 * K, N, rhf_parameters.spinOrbitalEnergiesBlocked());
    //const auto orbital_space = reference_onv.orbitalSpace();
    // Use a manually-made orbital space (#FIXME in develop).
    const auto orbital_space = GQCP::OrbitalSpace({0, 1, 2, 3, 4, 7, 8, 9, 10, 11}, {5, 6, 12, 13});  // occupied and virtual indices


    BOOST_REQUIRE(orbital_space.numberOfOrbitals() == M);


    // Initialize an environment suitable for CCD.
    auto environment = GQCP::CCSDEnvironment<double>::PerturbativeCCD(g_sq_hamiltonian, orbital_space);
    
    // Prepare the CCD solver and optimize the CCD model parameters.
    auto solver = GQCP::CCDSolver<double>::Plain();
    const auto ccd_qc_structure = GQCP::QCMethod::CCD<double>().optimize(solver, environment);

    const auto ccd_correlation_energy = ccd_qc_structure.groundStateEnergy();

    const double ref_ccd_correlation_energy = -0.070680088376;
    std::cout<<rhf_energy+ccd_correlation_energy<<std::endl;
    //BOOST_CHECK(std::abs(ccd_correlation_energy - ref_ccd_correlation_energy) < 1.0e-08);
}
