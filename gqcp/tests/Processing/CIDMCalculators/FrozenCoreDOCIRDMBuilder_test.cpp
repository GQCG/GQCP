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

#define BOOST_TEST_MODULE "FrozenCoreDOCIRDM_test"

#include <boost/test/unit_test.hpp>

#include "DensityMatrix/CIDMCalculators/FrozenCoreDOCIRDMBuilder.hpp"
#include "DensityMatrix/CIDMCalculators/SeniorityZeroDMCalculator.hpp"
#include "DensityMatrix/CIDMCalculators/SpinResolvedSelectedDMCalculator.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/CI/HamiltonianBuilder/FrozenCoreDOCI.hpp"


BOOST_AUTO_TEST_CASE(FrozenCoreDOCI_one_DMs) {

    size_t K = 5;
    GQCP::Molecule H5 = GQCP::Molecule::HChain(K, 1.1);
    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spinor_basis {H5, "STO-3G"};
    auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spinor_basis, H5);  // in an AO basis

    GQCP::SpinUnresolvedFrozenONVBasis fock_space {K, 3, 2};
    GQCP::SpinResolvedSelectedONVBasis selected_fock_space {fock_space};
    GQCP::FrozenCoreDOCI doci {fock_space};

    // Specify solver options and solve the eigenvalue problem
    // Solve the dense DOCI eigenvalue problem
    GQCP::CISolver ci_solver {doci, sq_hamiltonian};
    GQCP::DenseSolverOptions solver_options;
    ci_solver.solve(solver_options);

    GQCP::VectorX<double> coef = ci_solver.eigenpair().eigenvector();

    // Get the frozen core DOCI and SelectedCI 1-DMs
    GQCP::SpinResolvedSelectedDMCalculator sci_rdm {selected_fock_space};
    GQCP::FrozenCoreDOCIRDMBuilder doci_rdm {fock_space};
    GQCP::SpinResolved1DM<double> one_DMs_s = sci_rdm.calculateSpinResolved1DM(coef);
    GQCP::SpinResolved1DM<double> one_DMs = doci_rdm.calculateSpinResolved1DM(coef);

    BOOST_CHECK(one_DMs_s.orbitalDensity().isApprox(one_DMs.orbitalDensity()));
    BOOST_CHECK(one_DMs_s.alpha().isApprox(one_DMs.alpha()));
    BOOST_CHECK(one_DMs_sbeta().isApprox(one_DMs.beta()));
}


BOOST_AUTO_TEST_CASE(FrozenCoreDOCI_two_DMs) {

    size_t K = 5;
    GQCP::Molecule H5 = GQCP::Molecule::HChain(K, 1.1);
    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spinor_basis {H5, "STO-3G"};
    auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spinor_basis, H5);  // in an AO basis

    GQCP::SpinUnresolvedFrozenONVBasis fock_space {K, 3, 2};
    GQCP::SpinResolvedSelectedONVBasis selected_fock_space {fock_space};
    GQCP::FrozenCoreDOCI doci {fock_space};

    // Specify solver options and solve the eigenvalue problem
    // Solve the dense DOCI eigenvalue problem
    GQCP::CISolver ci_solver {doci, sq_hamiltonian};
    GQCP::DenseSolverOptions solver_options;
    ci_solver.solve(solver_options);

    GQCP::VectorX<double> coef = ci_solver.eigenpair().eigenvector();

    // Get the frozen core DOCI and SelectedCI 2-DMs
    GQCP::SpinResolvedSelectedDMCalculator sci_rdm {selected_fock_space};
    GQCP::FrozenCoreDOCIRDMBuilder doci_rdm {fock_space};
    GQCP::SpinResolved2DM<double> two_DMs_s = sci_rdm.calculateSpinResolved2DM(coef);
    GQCP::SpinResolved2DM<double> two_DMs = doci_rdm.calculateSpinResolved2DM(coef);

    BOOST_CHECK(two_DMs_s.alphaAlpha().isApprox(two_DMs.alphaAlpha(), 1.0e-06));
    BOOST_CHECK(two_DMs_s.alphaBeta().isApprox(two_DMs.alphaBeta(), 1.0e-06));
    BOOST_CHECK(two_DMs_s.betaAlpha().isApprox(two_DMs.betaAlpha(), 1.0e-06));
    BOOST_CHECK(two_DMs_s.betaBeta().isApprox(two_DMs.betaBeta(), 1.0e-06));
    BOOST_CHECK(two_DMs_s.orbitalDensity().isApprox(two_DMs.orbitalDensity(), 1.0e-06));
}
