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

#define BOOST_TEST_MODULE "SpinUnresolvedSelectedONVBasis"

#include <boost/test/unit_test.hpp>

#include "Basis/SpinorBasis/GSpinorBasis.hpp"
#include "Molecule/Molecule.hpp"
#include "ONVBasis/SpinUnresolvedSelectedONVBasis.hpp"


/**
 *  Check if the evaluation of a `GSQHamiltonian` in a `SpinUnresolvedSelectedONVBasis` works as expected, by diagonalizing a real and a complex Hamiltonian that differ by a complex unitary rotation.
 * 
 *  The eigenvalues found by the respective diagonalizations should be equal.
 */
BOOST_AUTO_TEST_CASE(GSQHamiltonian_rotation) {

    // Create the molecular Hamiltonian in the real Löwdin basis.
    const auto molecule = GQCP::Molecule::HChain(4, 1.0);
    GQCP::GSpinorBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    const auto M = spinor_basis.numberOfSpinors();
    spinor_basis.lowdinOrthonormalize();
    const auto hamiltonian = spinor_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));

    // Set up the full spin-unresolved selected ONV basis.
    const GQCP::SpinUnresolvedONVBasis onv_basis {M, molecule.numberOfElectrons()};
    const GQCP::SpinUnresolvedSelectedONVBasis selected_onv_basis {onv_basis};

    // Diagonalize the real Hamiltonian matrix representation.
    const auto H = selected_onv_basis.evaluateOperatorDense(hamiltonian);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver {H};


    // Create the molecular Hamiltonian in the Löwdin basis, and use a random complex unitary rotation.
    GQCP::GSpinorBasis<GQCP::complex, GQCP::GTOShell> spinor_basis_cd {molecule, "STO-3G"};
    spinor_basis_cd.lowdinOrthonormalize();
    const auto U = GQCP::GTransformation<GQCP::complex>::RandomUnitary(M);
    spinor_basis_cd.transform(U);
    const auto hamiltonian_cd = spinor_basis_cd.quantize(GQCP::FQMolecularHamiltonian(molecule));

    // Diagonalize the complex Hamiltonian matrix representation.
    const auto H_cd = selected_onv_basis.evaluateOperatorDense(hamiltonian_cd);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigensolver_cd {H_cd};


    BOOST_CHECK(eigensolver.eigenvalues().isApprox(eigensolver_cd.eigenvalues(), 1.0e-12));
}
