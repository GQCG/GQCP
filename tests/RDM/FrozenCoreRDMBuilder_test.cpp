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
#define BOOST_TEST_MODULE "Selected_RDM_test"


#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>

#include "RDM/RDMCalculator.hpp"
#include "RDM/SelectedRDMBuilder.hpp"

#include "RDM/DOCIRDMBuilder.hpp"
#include "RDM/FCIRDMBuilder.hpp"

#include "CISolver/CISolver.hpp"
#include "HamiltonianBuilder/FrozenCoreFCI.hpp"
#include "HamiltonianBuilder/DOCI.hpp"
#include "HamiltonianParameters/HamiltonianParameters.hpp"

#include "utilities/linalg.hpp"


BOOST_AUTO_TEST_CASE ( wakaka ) {

    size_t K = 6;
    GQCP::Molecule H5 = GQCP::Molecule::HChain(K, 1.1);
    auto ham_par = GQCP::HamiltonianParameters::Molecular(H5, "STO-3G");

    GQCP::FrozenProductFockSpace fock_space (K, 4, 4, 3);  // dim = 16
    GQCP::SelectedFockSpace sfock (fock_space);
    GQCP::FrozenCoreFCI fci (fock_space);

    // Specify solver options and solve the eigenvalue problem
    // Solve the dense FCI eigenvalue problem
    GQCP::CISolver ci_solver (fci, ham_par);
    GQCP::DenseSolverOptions solver_options;
    ci_solver.solve(solver_options);

    Eigen::VectorXd coef = ci_solver.get_eigenpair().get_eigenvector();

    // Get the 1-RDM from FCI
    GQCP::SelectedRDMBuilder sci_rdm(sfock);
    GQCP::OneRDMs one_rdms = sci_rdm.calculate1RDMs(coef);
    GQCP::TwoRDMs two_rdms = sci_rdm.calculate2RDMs(coef);




    std::cout<<one_rdms.one_rdm.get_matrix_representation()<<std::endl<<"------------------"<<std::endl;

    auto xy = two_rdms.two_rdm.get_matrix_representation();

    for (int i = 0; i < 3; i++) {
        for (int j = 3; j < K; j++) {
            for (int l = 3; l < K; l++) {
                for (int k = 0; k < 3; k++) {
                    std::cout<<i<<" "<<j<<" "<<l<<" "<<k<<" : "<<xy(i,j,l,k)<<std::endl;
                }
            }
        }
    }

    BOOST_CHECK(true);
}


