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

#define BOOST_TEST_MODULE "UHFSCFSolver"

#include <boost/test/unit_test.hpp>

#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/HF/UHF/UHFSCFSolver.hpp"


/**
 *  The following tests check if the default implementations for the UHF SCF solvers give the correct result.
 */


/**
 *  Check if our plain UHF SCF solver finds results (energy, orbital energies and coefficient matrix) that are equal to those from our RHF SCF algorithm.
 */
BOOST_AUTO_TEST_CASE(h2o_sto3g_plain) {

    // List the reference data
    const double ref_total_energy = -74.942080055631;

    GQCP::VectorX<double> ref_orbital_energies {7};  // the STO-3G basisset has 7 basis functions for water
    ref_orbital_energies << -20.26289322, -1.20969863, -0.54796582, -0.43652631, -0.38758791, 0.47762043, 0.5881361;

    GQCP::TransformationMatrix<double> ref_C_matrix {7};
    // clang-format off
    ref_C_matrix << -9.94434594e-01, -2.39158997e-01,  3.61117086e-17, -9.36837259e-02,  3.73303682e-31, -1.11639152e-01, -9.04958229e-17,
                    -2.40970260e-02,  8.85736467e-01, -1.62817254e-16,  4.79589270e-01, -1.93821120e-30,  6.69575233e-01,  5.16088339e-16,
                     1.59542752e-18,  5.29309704e-17, -6.07288675e-01, -1.49717339e-16,  8.94470461e-17, -8.85143477e-16,  9.19231270e-01,
                    -3.16155527e-03,  8.58957413e-02,  2.89059171e-16, -7.47426286e-01,  2.81871324e-30,  7.38494291e-01,  6.90314422e-16,
                     6.65079968e-35,  1.16150362e-32, -2.22044605e-16, -4.06685146e-30, -1.00000000e+00, -1.78495825e-31,  2.22044605e-16,
                     4.59373756e-03,  1.44038811e-01, -4.52995183e-01, -3.29475784e-01,  2.16823939e-16, -7.09847234e-01, -7.32462496e-01,
                     4.59373756e-03,  1.44038811e-01,  4.52995183e-01, -3.29475784e-01, -2.16823939e-16, -7.09847234e-01,  7.32462496e-01;
    // clang-format on
    const GQCP::UTransformationComponent<double> ref_C {ref_C_matrix};

    // Do our own UHF calculation.
    const auto water = GQCP::Molecule::ReadXYZ("data/h2o.xyz");
    const auto N_alpha = water.numberOfElectronPairs();
    const auto N_beta = water.numberOfElectronPairs();

    const GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spinor_basis {water, "STO-3G"};
    const auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spinor_basis, water);  // in an AO basis

    auto uhf_environment = GQCP::UHFSCFEnvironment<double>::WithCoreGuess(N_alpha, N_beta, sq_hamiltonian, spinor_basis.overlap().parameters());
    auto plain_uhf_scf_solver = GQCP::UHFSCFSolver<double>::Plain();
    plain_uhf_scf_solver.perform(uhf_environment);


    // Check the calculated results with the reference
    const double total_energy = uhf_environment.electronic_energies.back() + GQCP::Operator::NuclearRepulsion(water).value();
    BOOST_CHECK(std::abs(total_energy - ref_total_energy) < 1.0e-06);

    BOOST_CHECK(ref_orbital_energies.areEqualEigenvaluesAs(uhf_environment.orbital_energies_alpha.back(), 1.0e-06));
    BOOST_CHECK(ref_orbital_energies.areEqualEigenvaluesAs(uhf_environment.orbital_energies_beta.back(), 1.0e-06));

    BOOST_CHECK(ref_C.hasEqualSetsOfEigenvectorsAs(uhf_environment.coefficient_matrices.back().alpha(), 1.0e-05));
    BOOST_CHECK(ref_C.hasEqualSetsOfEigenvectorsAs(uhf_environment.coefficient_matrices.back().beta(), 1.0e-05));
}


/**
 *  Check if our DIIS UHF SCF solver finds results (energy, orbital energies and coefficient matrix) that are equal to those from our RHF SCF algorithm.
 */
BOOST_AUTO_TEST_CASE(h2o_sto3g_diis) {

    // List the reference data
    const double ref_total_energy = -74.942080055631;

    GQCP::VectorX<double> ref_orbital_energies {7};  // the STO-3G basisset has 7 basis functions for water
    ref_orbital_energies << -20.26289322, -1.20969863, -0.54796582, -0.43652631, -0.38758791, 0.47762043, 0.5881361;

    GQCP::TransformationMatrix<double> ref_C_matrix {7};
    // clang-format off
    ref_C_matrix << -9.94434594e-01, -2.39158997e-01,  3.61117086e-17, -9.36837259e-02,  3.73303682e-31, -1.11639152e-01, -9.04958229e-17,
                    -2.40970260e-02,  8.85736467e-01, -1.62817254e-16,  4.79589270e-01, -1.93821120e-30,  6.69575233e-01,  5.16088339e-16,
                     1.59542752e-18,  5.29309704e-17, -6.07288675e-01, -1.49717339e-16,  8.94470461e-17, -8.85143477e-16,  9.19231270e-01,
                    -3.16155527e-03,  8.58957413e-02,  2.89059171e-16, -7.47426286e-01,  2.81871324e-30,  7.38494291e-01,  6.90314422e-16,
                     6.65079968e-35,  1.16150362e-32, -2.22044605e-16, -4.06685146e-30, -1.00000000e+00, -1.78495825e-31,  2.22044605e-16,
                     4.59373756e-03,  1.44038811e-01, -4.52995183e-01, -3.29475784e-01,  2.16823939e-16, -7.09847234e-01, -7.32462496e-01,
                     4.59373756e-03,  1.44038811e-01,  4.52995183e-01, -3.29475784e-01, -2.16823939e-16, -7.09847234e-01,  7.32462496e-01;
    // clang-format on
    const GQCP::UTransformationComponent<double> ref_C {ref_C_matrix};

    // Do our own UHF calculation.
    const auto water = GQCP::Molecule::ReadXYZ("data/h2o.xyz");
    const auto N_alpha = water.numberOfElectronPairs();
    const auto N_beta = water.numberOfElectronPairs();

    const GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spinor_basis {water, "STO-3G"};
    const auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spinor_basis, water);  // in an AO basis

    auto uhf_environment = GQCP::UHFSCFEnvironment<double>::WithCoreGuess(N_alpha, N_beta, sq_hamiltonian, spinor_basis.overlap().parameters());
    auto diis_uhf_scf_solver = GQCP::UHFSCFSolver<double>::DIIS();
    diis_uhf_scf_solver.perform(uhf_environment);


    // Check the calculated results with the reference
    const double total_energy = uhf_environment.electronic_energies.back() + GQCP::Operator::NuclearRepulsion(water).value();
    BOOST_CHECK(std::abs(total_energy - ref_total_energy) < 1.0e-06);

    BOOST_CHECK(ref_orbital_energies.areEqualEigenvaluesAs(uhf_environment.orbital_energies_alpha.back(), 1.0e-06));
    BOOST_CHECK(ref_orbital_energies.areEqualEigenvaluesAs(uhf_environment.orbital_energies_beta.back(), 1.0e-06));

    BOOST_CHECK(ref_C.hasEqualSetsOfEigenvectorsAs(uhf_environment.coefficient_matrices.back().alpha(), 1.0e-05));
    BOOST_CHECK(ref_C.hasEqualSetsOfEigenvectorsAs(uhf_environment.coefficient_matrices.back().beta(), 1.0e-05));
}
