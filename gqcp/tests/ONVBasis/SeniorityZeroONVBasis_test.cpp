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

#define BOOST_TEST_MODULE "SeniorityZeroONVBasis"

#include <boost/test/unit_test.hpp>

#include "ONVBasis/SeniorityZeroONVBasis.hpp"
#include "ONVBasis/SpinResolvedSelectedONVBasis.hpp"


/**
 *  Check if the specific implementation of an operator's diagonal representation for seniority-zero ONV bases is equal to the generic selected ONV basis.
 */
BOOST_AUTO_TEST_CASE(evaluateOperatorDiagonal) {

    // Construct a molecular Hamiltonian in an orthonormal spinor basis.
    const auto molecule = GQCP::Molecule::HChain(3, 1.0);  // a H3-chain
    const auto N_P = molecule.numberOfElectronPairs();

    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    const auto K = spinor_basis.numberOfSpatialOrbitals();
    spinor_basis.lowdinOrthonormalize();

    const auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spinor_basis, molecule);

    // Set up the specific and generic ONV bases.
    const GQCP::SeniorityZeroONVBasis sz_onv_basis {K, N_P};
    const GQCP::SpinResolvedSelectedONVBasis selected_onv_basis {sz_onv_basis};


    // Check the if the evaluations are correct.
    const auto& h = sq_hamiltonian.core();
    const auto& g = sq_hamiltonian.twoElectron();

    BOOST_CHECK(sz_onv_basis.evaluateOperatorDiagonal(h).isApprox(selected_onv_basis.evaluateOperatorDiagonal(h), 1.0e-08));
    BOOST_CHECK(sz_onv_basis.evaluateOperatorDiagonal(g).isApprox(selected_onv_basis.evaluateOperatorDiagonal(g), 1.0e-08));
    BOOST_CHECK(sz_onv_basis.evaluateOperatorDiagonal(sq_hamiltonian).isApprox(selected_onv_basis.evaluateOperatorDiagonal(sq_hamiltonian), 1.0e-08));
}
