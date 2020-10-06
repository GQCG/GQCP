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

#define BOOST_TEST_MODULE "SpinResolvedFrozenDMCalculator_test"

#include <boost/test/unit_test.hpp>

#include "DensityMatrix/CIDMCalculators/SpinResolvedFrozenDMCalculator.hpp"
#include "DensityMatrix/CIDMCalculators/SpinResolvedSelectedDMCalculator.hpp"
#include "Mathematical/Optimization/Eigenproblem/EigenproblemSolver.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/CI/CI.hpp"
#include "QCMethod/CI/CIEnvironment.hpp"


/**
 *  Check if the 1- and 2-DMs for a frozen core spin-resolved ONV basis are equal to the 'selected' case.
 *  The system of interest is a linear chain of 5 H atoms, 1.1 bohr apart, using an STO-3G basisset.
 */
BOOST_AUTO_TEST_CASE(FrozenCoreFCI_one_DMs) {

    // Set up the molecular Hamiltonian for H5//STO-3G in the Löwdin basis.
    const GQCP::Molecule molecule = GQCP::Molecule::HChain(5, 1.1);
    const auto N_alpha = 3;
    const auto N_beta = 2;

    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    const auto K = spinor_basis.numberOfSpatialOrbitals();
    spinor_basis.lowdinOrthonormalize();

    const auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, molecule);


    // Do a dense frozen core CI calculation for 2 frozen orbitals.
    const GQCP::SpinResolvedFrozenONVBasis onv_basis {K, N_alpha, N_beta, 2};

    auto environment = GQCP::CIEnvironment::Dense(sq_hamiltonian, onv_basis);
    auto solver = GQCP::EigenproblemSolver::Dense();

    const auto linear_expansion = GQCP::QCMethod::CI<GQCP::SpinResolvedFrozenONVBasis>(onv_basis).optimize(solver, environment).groundStateParameters();


    // Calculate the 1-DMs using specialized spin-resolved and 'selected' routines, and check if they are equal.
    const GQCP::SpinResolvedFrozenDMCalculator spin_resolved_dm_calculator {onv_basis};
    const auto one_DMs_specialized = spin_resolved_dm_calculator.calculateSpinResolved1DM(linear_expansion.coefficients());

    const GQCP::SpinResolvedSelectedONVBasis selected_onv_basis {onv_basis};
    const GQCP::SpinResolvedSelectedDMCalculator selected_dm_calculator {selected_onv_basis};
    const auto one_DMs_selected = selected_dm_calculator.calculateSpinResolved1DM(linear_expansion.coefficients());

    BOOST_CHECK(one_DMs_specialized.spinSummed().isApprox(one_DMs_selected.spinSummed(), 1.0e-12));
    BOOST_CHECK(one_DMs_specialized.alpha().isApprox(one_DMs_selected.alpha(), 1.0e-12));
    BOOST_CHECK(one_DMs_specialized.beta().isApprox(one_DMs_selected.beta(), 1.0e-12));


    // Calculate the 2-DMs using specialized spin-resolved and 'selected' routines, and check if they are equal.
    const auto two_DMs_specialized = spin_resolved_dm_calculator.calculateSpinResolved2DM(linear_expansion.coefficients());
    const auto two_DMs_selected = selected_dm_calculator.calculateSpinResolved2DM(linear_expansion.coefficients());

    BOOST_CHECK(two_DMs_specialized.alphaAlpha().isApprox(two_DMs_selected.alphaAlpha(), 1.0e-12));
    BOOST_CHECK(two_DMs_specialized.alphaBeta().isApprox(two_DMs_selected.alphaBeta(), 1.0e-12));
    BOOST_CHECK(two_DMs_specialized.betaAlpha().isApprox(two_DMs_selected.betaAlpha(), 1.0e-12));
    BOOST_CHECK(two_DMs_specialized.betaBeta().isApprox(two_DMs_selected.betaBeta(), 1.0e-12));
    BOOST_CHECK(two_DMs_specialized.spinSummed().isApprox(two_DMs_selected.spinSummed(), 1.0e-12));
}
