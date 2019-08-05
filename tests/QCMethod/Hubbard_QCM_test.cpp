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
#define BOOST_TEST_MODULE "HUBBARD_QCM_TEST"
#include <boost/test/unit_test.hpp>

#include "QCMethod/Hubbard.hpp"
#include "HamiltonianBuilder/Hubbard.hpp"
#include "CISolver/CISolver.hpp"
#include "RDM/RDMCalculator.hpp"


BOOST_AUTO_TEST_CASE ( Hubbard_QCM ) {

    // Check if the Hubbard module and the Hubbard QCM produce the same results
    // retrieve a suitable csline
    std::string csline = "-0.999984,-0.736924,0.511211,-0.082700,0.065534,-0.562082,-0.905911,0.357729,0.358593,0.869386";
    size_t K = 4;
    size_t nos = 2;  // number of states
    auto H = GQCP::HoppingMatrix::FromCSLine(csline);
    auto mol_ham_par = GQCP::HamiltonianParameters<double>::Hubbard(H);
    size_t N = 2;

    GQCP::QCMethod::Hubbard hubbard_qcm(csline, nos, N, N);
    hubbard_qcm.solve();
    auto const& hubbard_qcm_energies = hubbard_qcm.energies();
    auto const& hubbard_qcm_1rdms = hubbard_qcm.oneRDMs();

    // Create the Hubbard module
    GQCP::ProductFockSpace fock_space (K, N, N);  // dim = 36
    GQCP::Hubbard hubbard (fock_space);
    GQCP::RDMCalculator rdm_calc (fock_space);

    // Solve via dense
    GQCP::CISolver hubbard_solver (hubbard, mol_ham_par);

    GQCP::DenseSolverOptions dense_solver_options;
    dense_solver_options.number_of_requested_eigenpairs = nos;
    hubbard_solver.solve(dense_solver_options);

    auto const& hubbard_eigenpairs = hubbard_solver.get_eigenpairs();

    for (size_t i = 0; i < nos; i++) {
            BOOST_CHECK(std::abs(hubbard_qcm_energies[i] - (hubbard_eigenpairs[i].get_eigenvalue())) < 1.0e-06);
            rdm_calc.set_coefficients(hubbard_eigenpairs[i].get_eigenvector());
            auto const& one_rdm = rdm_calc.calculate1RDMs().one_rdm;
            BOOST_CHECK(one_rdm.isApprox(hubbard_qcm_1rdms[i], 1.0e-06));
    }
}

