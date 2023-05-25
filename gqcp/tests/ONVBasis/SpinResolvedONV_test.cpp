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

#define BOOST_TEST_MODULE "SpinResolvedONV"

#include <boost/test/unit_test.hpp>

#include "Basis/SpinorBasis/RSpinOrbitalBasis.hpp"
#include "Basis/SpinorBasis/USpinOrbitalBasis.hpp"
#include "ONVBasis/SpinResolvedONV.hpp"
#include "QCMethod/HF/RHF/DiagonalRHFFockMatrixObjective.hpp"
#include "QCMethod/HF/RHF/RHF.hpp"
#include "QCMethod/HF/RHF/RHFSCFSolver.hpp"
#include "QCModel/HF/UHF.hpp"


/**
 *  Check if the creation of the 'RHF' ONV works as expected.
 */
BOOST_AUTO_TEST_CASE(RHF) {

    // For K=5 alpha-spinorbitals, K=5 beta-spinorbitals, and N_P=3 electron electron pairs, the 'RHF' ONV should be ("00111" = 7) ("00111" = 7)
    const size_t K = 5;
    const size_t N_P = 3;
    const GQCP::SpinUnresolvedONV alpha_onv {K, N_P, 7};

    const GQCP::SpinResolvedONV reference {alpha_onv, alpha_onv};
    BOOST_CHECK(reference == GQCP::SpinResolvedONV::RHF(K, N_P));
}


/**
 *  Check if the creation of the 'UHF' ONV works as expected.
 */
BOOST_AUTO_TEST_CASE(UHF) {

    // For K=5 alpha-spinorbitals, K=5 beta-spinorbitals, and N_alpha=3 alpha-electrons and N_beta=2 beta-electrons, the 'UHF' ONV should be ("00111" = 7) ("00011" = 3)
    const size_t K = 5;
    const size_t N_alpha = 3;
    const size_t N_beta = 2;
    const GQCP::SpinUnresolvedONV alpha_onv {K, N_alpha, 7};
    const GQCP::SpinUnresolvedONV beta_onv {K, N_beta, 3};

    const GQCP::SpinResolvedONV reference {alpha_onv, beta_onv};
    BOOST_CHECK(reference == GQCP::SpinResolvedONV::UHF(K, N_alpha, N_beta));
}


/**
 *  Check if the creation of a spin-resolved ONV from a textual/string representation works as expected.
 */
BOOST_AUTO_TEST_CASE(FromString) {

    const GQCP::SpinUnresolvedONV onv_alpha_ref {5, 3, 21};  // "10101" (21)
    const GQCP::SpinUnresolvedONV onv_beta_ref {5, 3, 22};   // "10110" (22)
    const GQCP::SpinResolvedONV onv_ref {onv_alpha_ref, onv_beta_ref};

    const auto onv = GQCP::SpinResolvedONV::FromString("10101", "10110");
    BOOST_CHECK(onv == onv_ref);
}


/**
 *  Check if the overlap between an RHF-related ONV and an UHF-related ONV works as expected.
 *
 *  We don't really have a reference implementation, but we can check if the overlaps are equal to 1 or 0 if we use RHF orbitals and UHF orbitals that have the same alpha- and beta-part.
 *
 *  The system under consideration is H2 with a STO-3G basisset.
 */
BOOST_AUTO_TEST_CASE(RHF_UHF_overlap) {

    // Create an RSpinOrbitalBasis with the canonical RHF spin-orbitals.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2.xyz");
    const auto N = molecule.numberOfElectrons();
    const auto N_P = molecule.numberOfElectronPairs();

    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> r_spinor_basis {molecule, "STO-3G"};
    const auto K = r_spinor_basis.numberOfSpatialOrbitals();
    const auto S = r_spinor_basis.overlap().parameters();

    auto sq_hamiltonian = r_spinor_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));  // in an AO basis

    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(N, sq_hamiltonian, S);
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective {sq_hamiltonian};

    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters();

    r_spinor_basis.transform(rhf_parameters.expansion());


    // Convert the RSpinOrbitalBasis into an USpinOrbitalBasis, yielding an unrestricted spin-orbital basis where C_alpha == C_beta.
    const auto u_spinor_basis = GQCP::USpinOrbitalBasis<double, GQCP::GTOShell>::FromRestricted(r_spinor_basis);


    // Check if the UHF determinant has overlap 1 with the corresponding RHF determinant, and overlap 0 with the excitations on top of the RHF reference.
    const auto& C_restricted = rhf_parameters.expansion();
    const auto C_unrestricted = GQCP::UTransformation<double>::FromRestricted(C_restricted);

    const auto uhf_determinant = GQCP::SpinResolvedONV::UHF(K, N_P, N_P);

    const auto rhf_onv_0101 = GQCP::SpinResolvedONV::FromString("01", "01");
    const auto rhf_onv_0110 = GQCP::SpinResolvedONV::FromString("01", "10");
    const auto rhf_onv_1001 = GQCP::SpinResolvedONV::FromString("10", "01");
    const auto rhf_onv_1010 = GQCP::SpinResolvedONV::FromString("10", "10");

    BOOST_CHECK(std::abs(uhf_determinant.calculateProjection(rhf_onv_0101, C_unrestricted, C_restricted, S) - 1.0) < 1.0e-12);
    BOOST_CHECK(std::abs(uhf_determinant.calculateProjection(rhf_onv_0110, C_unrestricted, C_restricted, S) - 0.0) < 1.0e-12);
    BOOST_CHECK(std::abs(uhf_determinant.calculateProjection(rhf_onv_1001, C_unrestricted, C_restricted, S) - 0.0) < 1.0e-12);
    BOOST_CHECK(std::abs(uhf_determinant.calculateProjection(rhf_onv_1010, C_unrestricted, C_restricted, S) - 0.0) < 1.0e-12);
}
