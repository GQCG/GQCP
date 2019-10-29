// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2019  the GQCG developers
// 
// GQCG-gqcp is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// GQCG-gqcp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-gqcp.  If not, see <http://www.gnu.org/licenses/>.
// 
#define BOOST_TEST_MODULE "SPBasis_test"

#include <boost/test/unit_test.hpp>

#include "Basis/SpinorBasis/RSpinorBasis.hpp"

#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "RHF/PlainRHFSCFSolver.hpp"


/**
 *  Check if using a Löwdin orthonormalization ensures an orthonormal single-particle basis
 */
BOOST_AUTO_TEST_CASE ( Lowdin_orthonormal ) {

    // Construct the initial single-particle basis (corresponding to the underlying GTOs)
    const auto h2 = GQCP::Molecule::ReadXYZ("data/h2.xyz");
    GQCP::RSpinorBasis<double, GQCP::GTOShell> sp_basis (h2, "STO-3G");
    BOOST_REQUIRE_EQUAL(sp_basis.numberOfOrbitals(), 2);


    // Löwdin-orthonormalize and check the result
    sp_basis.lowdinOrthonormalize();
    BOOST_CHECK(sp_basis.isOrthonormal());
}


/**
 *  Check if the Löwdin-orthonormalization matrix depends on the current orbitals: the Löwdin basis isn't reached when T=S_AO^{-1/2}, but when T=S_current^{-1/2}
 */
BOOST_AUTO_TEST_CASE ( lowdinOrthonormalizatioMatrix ) {

    // Construct the initial single-particle basis (corresponding to the underlying GTOs) and calculate the corresponding Löwdin orthonormalization matrix
    const auto h2 = GQCP::Molecule::ReadXYZ("data/h2.xyz");
    GQCP::RSpinorBasis<double, GQCP::GTOShell> sp_basis (h2, "STO-3G");
    const auto T_lowdin_1 = sp_basis.lowdinOrthonormalizationMatrix();


    // Transform the single-particle basis and re-calculate the Löwdin orthonormalization matrix and check the result
    const auto K = sp_basis.numberOfBasisFunctions();
    sp_basis.transform(GQCP::TransformationMatrix<double>::Random(K, K));
    const auto T_lowdin_2 = sp_basis.lowdinOrthonormalizationMatrix();

    BOOST_CHECK(!T_lowdin_1.isApprox(T_lowdin_2, 1.0e-08));  // the two Löwdin transformation matrices should not be equal
}


/**
 *  Check if the orbitals in an AO basis are not orthonormal, but after a transformation to the canonical RHF orbitals, they are
 */
BOOST_AUTO_TEST_CASE ( isOrthonormal ) {

    // The orbitals in an AO basis are not orthonormal
    const auto h2o = GQCP::Molecule::ReadXYZ("data/h2o.xyz");
    GQCP::RSpinorBasis<double, GQCP::GTOShell> sp_basis (h2o, "STO-3G");
    BOOST_CHECK(!sp_basis.isOrthonormal());


    // The orbitals in the RHF basis should be orthonormal
    const auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(sp_basis, h2o);  // in the AO basis
    GQCP::PlainRHFSCFSolver plain_scf_solver (sq_hamiltonian, sp_basis, h2o);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();
    sp_basis.transform(rhf.get_C());

    BOOST_CHECK(sp_basis.isOrthonormal());
}
