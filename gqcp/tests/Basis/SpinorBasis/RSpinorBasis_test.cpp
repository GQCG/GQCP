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
#define BOOST_TEST_MODULE "SimpleSpinorBasis_test"

#include <boost/test/unit_test.hpp>

#include "Basis/SpinorBasis/RSpinorBasis.hpp"

#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/RHF/PlainRHFSCFSolver.hpp"


/**
 *  Check if the orbitals in an AO basis are not orthonormal, but after a transformation to the canonical RHF orbitals, they are
 */
BOOST_AUTO_TEST_CASE ( RHF_orbitals_are_orthonormal ) {

    // The orbitals in an AO basis are not orthonormal
    const auto h2o = GQCP::Molecule::ReadXYZ("data/h2o.xyz");
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (h2o, "STO-3G");
    BOOST_CHECK(!spinor_basis.isOrthonormal());


    // The orbitals in the RHF basis should be orthonormal
    const auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, h2o);  // in the AO basis
    GQCP::PlainRHFSCFSolver plain_scf_solver (sq_hamiltonian, spinor_basis, h2o);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();
    spinor_basis.transform(rhf.get_C());

    BOOST_CHECK(spinor_basis.isOrthonormal());
}
