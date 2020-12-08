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

#define BOOST_TEST_MODULE "Simple2DM"

#include <boost/test/unit_test.hpp>

#include "DensityMatrix/G2DM.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/HF/RHF/DiagonalRHFFockMatrixObjective.hpp"
#include "QCMethod/HF/RHF/RHF.hpp"
#include "QCMethod/HF/RHF/RHFSCFSolver.hpp"


/*
 *  MARK: Helper functions
 */

/**
 *  @return a toy 2-DM where
 *      d(i,j,k,l) = l + 2k + 4j + 8i
 */
GQCP::G2DM<double> calculateToy2DMTensor() {

    GQCP::G2DM<double> d {2};

    for (size_t i = 0; i < 2; i++) {
        for (size_t j = 0; j < 2; j++) {
            for (size_t k = 0; k < 2; k++) {
                for (size_t l = 0; l < 2; l++) {
                    auto i_ = static_cast<double>(i);
                    auto j_ = static_cast<double>(j);
                    auto k_ = static_cast<double>(k);
                    auto l_ = static_cast<double>(l);

                    d(i, j, k, l) = l_ + 2 * k_ + 4 * j_ + 8 * i_;
                }
            }
        }
    }

    return d;
};


/*
 *  MARK: Unit tests
 */


/**
 *  Check if the 2-DM `trace` method is correctly implemented, from a manual calculation.
 */
BOOST_AUTO_TEST_CASE(trace) {

    const auto d = calculateToy2DMTensor();

    BOOST_CHECK(std::abs(d.trace() - 30.0) < 1.0e-12);
}


/**
 *  Check if the 2-DM `reduce` method is correctly implemented, from a manual calculation.
 */
BOOST_AUTO_TEST_CASE(reduce) {

    const auto d = calculateToy2DMTensor();

    // Set up the reference result.
    GQCP::G1DM<double> D_ref = GQCP::G1DM<double>::Zero(2);

    // clang-format off
    D_ref <<  3, 11,
             19, 27;
    // clang-format on

    BOOST_CHECK(D_ref.isApprox(d.reduce(), 1.0e-12));
}


/**
 *  Test if the expectation value of a two-electron operator in different orbital bases is the same.
 */
BOOST_AUTO_TEST_CASE(two_electron_operator_expectation_value_different_orbital_bases) {

    // Prepare the molecular Hamiltonian in the AO basis, in order to proceed with an RHF SCF calculation.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/ch4_crawdad.xyz");
    const auto N = molecule.numberOfElectrons();

    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "STO-3G"};
    const auto S = spin_orbital_basis.overlap();

    const auto hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spin_orbital_basis, molecule);  // In the AO basis.
    const auto K = hamiltonian.numberOfOrbitals();

    // Do the RHF SCF calculation to retrieve the RHF MOs.
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(molecule.numberOfElectrons(), hamiltonian, S.parameters());
    auto diis_rhf_scf_solver = GQCP::RHFSCFSolver<double>::DIIS();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective {hamiltonian};

    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, diis_rhf_scf_solver, rhf_environment).groundStateParameters();


    // Prepare three two-electron operators in different orbital basis.
    const auto& g_AO = hamiltonian.twoElectron();
    const auto g_MO = g_AO.transformed(rhf_parameters.expansion());

    const auto T_random = GQCP::RTransformation<double>::Random(K);
    const auto g_random = g_AO.transformed(T_random);

    // Prepare three density matrices in the corresponding orbital bases.
    const auto d_AO = rhf_parameters.calculateOrthonormalBasis2DM();
    const auto d_MO = d_AO.transformed(rhf_parameters.expansion());
    const auto d_random = d_AO.transformed(T_random);


    // Check if the expectation values match.
    const double exp_val_AO = g_AO.calculateExpectationValue(d_AO);
    const double exp_val_MO = g_MO.calculateExpectationValue(d_MO);
    const double exp_val_random = g_random.calculateExpectationValue(d_random);

    BOOST_CHECK(std::abs(exp_val_AO - exp_val_MO) < 1.0e-11);
    BOOST_CHECK(std::abs(exp_val_AO - exp_val_random) < 1.0e-11);
}
