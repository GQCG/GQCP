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

#define BOOST_TEST_MODULE "SeniorityZeroDMCalculator_test"

#include <boost/test/unit_test.hpp>

#include "Mathematical/Optimization/Eigenproblem/EigenproblemSolver.hpp"
#include "ONVBasis/SeniorityZeroONVBasis.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "Processing/DensityMatrices/CIDMCalculators/SeniorityZeroDMCalculator.hpp"
#include "Processing/DensityMatrices/CIDMCalculators/SpinResolvedSelectedDMCalculator.hpp"
#include "QCMethod/CI/CI.hpp"
#include "QCMethod/CI/CIEnvironment.hpp"


/**
 *  Check if the 1- and 2-DMs for a seniority-zero ONV basis are equal to the 'selected' case.
 *  The system of interested is H2O//STO-3G, with 7 spatial orbitals and a Fock space dimension of 441.
 */
BOOST_AUTO_TEST_CASE(specialized_vs_selected_DMs) {

    // Set up the molecular Hamiltonian in a Löwdin-orthonormalized spinor basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o_Psi4_GAMESS.xyz");
    const auto N_P = molecule.numberOfElectrons() / 2;  // number of electron pairs

    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    const auto K = spinor_basis.numberOfSpatialOrbitals();
    spinor_basis.lowdinOrthonormalize();

    const auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, molecule);


    // Do a dense DOCI calculation.
    const GQCP::SeniorityZeroONVBasis onv_basis {K, N_P};

    auto environment = GQCP::CIEnvironment::Dense(sq_hamiltonian, onv_basis);
    auto solver = GQCP::EigenproblemSolver::Dense();

    const auto linear_expansion = GQCP::QCMethod::CI<GQCP::SeniorityZeroONVBasis>(onv_basis).optimize(solver, environment).groundStateParameters();


    // Calculate the 1-DMs using specialized spin-resolved and 'selected' routines, and check if they are equal.
    const GQCP::SeniorityZeroDMCalculator seniority_zero_dm_calculator {onv_basis};
    const auto one_DMs_specialized = seniority_zero_dm_calculator.calculateSpinResolved1DM(linear_expansion.coefficients());

    const GQCP::SpinResolvedSelectedONVBasis selected_onv_basis {onv_basis};
    const GQCP::SpinResolvedSelectedDMCalculator selected_dm_calculator {selected_onv_basis};
    const auto one_DMs_selected = selected_dm_calculator.calculateSpinResolved1DM(linear_expansion.coefficients());

    BOOST_CHECK(one_DMs_specialized.spinSummed().isApprox(one_DMs_selected.spinSummed(), 1.0e-12));
    BOOST_CHECK(one_DMs_specialized.alpha().isApprox(one_DMs_selected.alpha(), 1.0e-12));
    BOOST_CHECK(one_DMs_specialized.beta().isApprox(one_DMs_selected.beta(), 1.0e-12));


    // Calculate the 2-DMs using specialized spin-resolved and 'selected' routines, and check if they are equal.
    const auto two_DMs_specialized = seniority_zero_dm_calculator.calculateSpinResolved2DM(linear_expansion.coefficients());
    const auto two_DMs_selected = selected_dm_calculator.calculateSpinResolved2DM(linear_expansion.coefficients());

    BOOST_CHECK(two_DMs_specialized.alphaAlpha().isApprox(two_DMs_selected.alphaAlpha(), 1.0e-12));
    BOOST_CHECK(two_DMs_specialized.alphaBeta().isApprox(two_DMs_selected.alphaBeta(), 1.0e-12));
    BOOST_CHECK(two_DMs_specialized.betaAlpha().isApprox(two_DMs_selected.betaAlpha(), 1.0e-12));
    BOOST_CHECK(two_DMs_specialized.betaBeta().isApprox(two_DMs_selected.betaBeta(), 1.0e-12));
    BOOST_CHECK(two_DMs_specialized.spinSummed().isApprox(two_DMs_selected.spinSummed(), 1.0e-12));
}
