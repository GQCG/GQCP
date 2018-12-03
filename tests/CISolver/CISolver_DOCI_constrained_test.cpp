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
#define BOOST_TEST_MODULE "ConstrainedDociSolver"


#include "CISolver/CISolver.hpp"
#include "RDM/RDMCalculator.hpp"
#include "HamiltonianBuilder/DOCI.hpp"
#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "RHF/PlainRHFSCFSolver.hpp"
#include "RHF/DIISRHFSCFSolver.hpp"
#include "properties/expectation_values.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


BOOST_AUTO_TEST_CASE ( CO_DOCI_constrained_dense ) {
    // Reference data from tmhuysens thesis
    Eigen::Matrix<double, 21 , 3> CO_data;
    //        lamda | C population | DOCI energy
    CO_data <<  -1.000000, 	5.80863087698, 	-111.281035626,
                -0.900000, 	5.81016803318, 	-111.282495093,
                -0.800000, 	5.8118219348 , 	-111.28389979 ,
                -0.700000, 	5.81362985899, 	-111.285254278,
                -0.600000, 	5.81563610261, 	-111.286556475,
                -0.500000, 	5.81789473016, 	-111.28779636 ,
                -0.400000, 	5.82047326851, 	-111.288953711,
                -0.300000, 	5.82345780552, 	-111.289994495,
                -0.200000, 	5.82696018148, 	-111.290865222,
                -0.100000, 	5.83112832986, 	-111.291484159,
                 0.000000, 	5.83616141374, 	-111.291727603,
                 0.100000, 	5.84233235124, 	-111.291408185,
                 0.200000, 	5.85002183736, 	-111.29024015 ,
                 0.300000, 	5.85977033411, 	-111.28778306 ,
                 0.400000, 	5.87235795874, 	-111.283349644,
                 0.500000, 	5.88892639439, 	-111.275854656,
                 0.600000, 	5.91115885823, 	-111.263570726,
                 0.700000, 	5.94152162485, 	-111.243754333,
                 0.800000, 	5.98350566877, 	-111.212151981,
                 0.900000, 	6.0416010444 , 	-111.162616424,
                 1.000000, 	6.12030270656, 	-111.087663866;

    // Create the molecular Hamiltonian parameters for CO
    GQCP::Molecule CO ("../tests/data/CO_mulliken.xyz");
    auto mol_ham_par = GQCP::HamiltonianParameters::Molecular(CO, "STO-3G");

    // Create a plain RHF SCF solver and solve the SCF equations
    GQCP::DIISRHFSCFSolver diis_scf_solver (mol_ham_par, CO);
    diis_scf_solver.solve();
    auto rhf = diis_scf_solver.get_solution();

    // Transform the ham_par
    mol_ham_par.transform(rhf.get_C());

    GQCP::FockSpace fock_space (mol_ham_par.get_K(), CO.get_N()/2);  // dim = 4

    // Create the DOCI module
    GQCP::DOCI doci (fock_space);
    GQCP::RDMCalculator rdm_calc (fock_space);

    numopt::eigenproblem::DenseSolverOptions dense_solver_options;
    GQCP::Vectoru ao_list = {0,1,2,3,4};

    for (int i = 0; i < 21; i++) {
        // Calculate the Mulliken operator
        auto mulliken_operator = mol_ham_par.calculateMullikenOperator(ao_list);

        // Contrain the original Hamiltonian parameters
        auto constrained_ham_par = mol_ham_par.constrain(mulliken_operator, CO_data(i, 0));

        GQCP::CISolver ci_solver(doci, constrained_ham_par);
        // Solve Dense
        ci_solver.solve(dense_solver_options);

        // Retrieve the eigenvalues
        auto doci_dense_eigenvalue = ci_solver.get_eigenpair().get_eigenvalue();
        auto doci_dense_eigenvector = ci_solver.get_eigenpair().get_eigenvector();

        auto one_rdms = rdm_calc.calculate1RDMs(doci_dense_eigenvector);

        // Retrieve the expectation value of the Mulliken operator (aka the population)
        double mulliken_population = GQCP::calculateExpectationValue(mulliken_operator, one_rdms.one_rdm);

        // Retrieve the total energy by adding the lambda times the expectation value of the constraining operator
        double total_energy = doci_dense_eigenvalue + CO_data(i, 0) * mulliken_population + CO.calculateInternuclearRepulsionEnergy();


        BOOST_CHECK(std::abs(total_energy - CO_data(i, 2)) < 1.0e-3);
        BOOST_CHECK(std::abs(mulliken_population - CO_data(i, 1)) < 1.0e-2);
    }
    BOOST_CHECK(true);
}
