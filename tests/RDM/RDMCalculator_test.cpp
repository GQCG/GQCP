// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2018  the GQCG developers
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
#define BOOST_TEST_MODULE "RDMCalculator_test"


#include "RDM/RDMCalculator.hpp"

#include "FockSpace/FockSpace.hpp"
#include "FockSpace/ProductFockSpace.hpp"
#include "CISolver/CISolver.hpp"
#include "HamiltonianBuilder/DOCI.hpp"
#include "HamiltonianParameters/HamiltonianParameters.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


BOOST_AUTO_TEST_CASE ( constructor ) {

    // Test polymorphic entry for RDM (from DOCIRDMBuilder test-case).

    // Get the 1-RDM from DOCI
    size_t N = 4;  // 4 electrons
    auto ham_par = GQCP::HamiltonianParameters::ReadFCIDUMP("../tests/data/lih_631g_caitlin.FCIDUMP");
    size_t K = ham_par.get_K();  // 16 SO

    // Abstract pointer to test RDM
    std::shared_ptr<GQCP::BaseFockSpace> fock_space_dy(new GQCP::FockSpace(K, N/2));  // dim = 120
    GQCP::FockSpace fock_space (K, N/2);  // dim = 120

    GQCP::DOCI doci (fock_space);

    // Specify solver options and solve the eigenvalue problem
    // Solve the dense DOCI eigenvalue problem
    GQCP::CISolver ci_solver (doci, ham_par);
    GQCP::DenseSolverOptions solver_options;
    ci_solver.solve(solver_options);

    Eigen::VectorXd coef = ci_solver.get_eigenpair().get_eigenvector();
    
    // Check if the DOCI 1-RDM has the proper trace.
    GQCP::RDMCalculator doci_rdm (*fock_space_dy);
    GQCP::OneRDMs one_rdms = doci_rdm.calculate1RDMs(coef);

    BOOST_CHECK(std::abs(one_rdms.one_rdm.trace() - N) < 1.0e-12);
}
