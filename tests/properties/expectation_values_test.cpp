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
#define BOOST_TEST_MODULE "expectation_values"

#include "properties/expectation_values.hpp"

#include "RDM/RDMCalculator.hpp"
#include "CISolver/CISolver.hpp"
#include "HamiltonianBuilder/DOCI.hpp"
#include "HamiltonianParameters/HamiltonianParameters_constructors.hpp"
#include "RHF/DIISRHFSCFSolver.hpp"
#include "RHF/PlainRHFSCFSolver.hpp"
#include "LibintCommunicator.hpp"
#include "units.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>


BOOST_AUTO_TEST_CASE ( one_electron_throw ) {

    GQCP::OneElectronOperator h (Eigen::MatrixXd::Zero(2, 2));
    GQCP::OneRDM D_valid (Eigen::MatrixXd::Zero(2, 2));
    GQCP::OneRDM D_invalid (Eigen::MatrixXd::Zero(3, 3));

    BOOST_CHECK_THROW(GQCP::calculateExpectationValue(h, D_invalid), std::invalid_argument);
    BOOST_CHECK_NO_THROW(GQCP::calculateExpectationValue(h, D_valid));
}


BOOST_AUTO_TEST_CASE ( two_electron_throw ) {

    Eigen::Tensor<double, 4> g_tensor (2, 2, 2, 2);
    g_tensor.setZero();
    GQCP::TwoElectronOperator g (g_tensor);

    Eigen::Tensor<double, 4> d_tensor_valid (2, 2, 2, 2);
    d_tensor_valid.setZero();
    Eigen::Tensor<double, 4> d_tensor_invalid (3, 3, 3, 3);
    d_tensor_valid.setZero();
    GQCP::TwoRDM d_valid (d_tensor_valid);
    GQCP::TwoRDM d_invalid (d_tensor_invalid);

    BOOST_CHECK_THROW(GQCP::calculateExpectationValue(g, d_invalid), std::invalid_argument);
    BOOST_CHECK_NO_THROW(GQCP::calculateExpectationValue(g, d_valid));
}

BOOST_AUTO_TEST_CASE ( mulliken_N2_STO_3G ) {

    // Check that the mulliken population of N2 is 14 (N)

    // Initialize the molecule and molecular Hamiltonian parameters for N2
    GQCP::Atom N_1 (7, 0.0, 0.0, 0.0);
    GQCP::Atom N_2 (7, 0.0, 0.0, GQCP::units::angstrom_to_bohr(1.134));  // from CCCBDB, STO-3G geometry
    std::vector<GQCP::Atom> atoms {N_1, N_2};
    GQCP::Molecule N2 (atoms);

    auto ao_basis = std::make_shared<GQCP::AOBasis>(N2, "STO-3G");
    auto ham_par = GQCP::constructMolecularHamiltonianParameters(ao_basis);

    size_t K = ao_basis->get_number_of_basis_functions();

    // We include all basis functions
    GQCP::Vectoru gto_list (K);
    for(size_t i = 0; i<K; i++){
        gto_list[i] = i;
    }

    GQCP::OneElectronOperator mulliken = ham_par.calculateMullikenOperator(gto_list);

    size_t N = N2.get_N();

    // Create a 1-RDM for N2
    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(K, K);
    GQCP::OneRDM one_rdm = GQCP::calculateRHF1RDM(K, N);

    double mulliken_population = GQCP::calculateExpectationValue(mulliken, one_rdm);
    BOOST_CHECK(std::abs(mulliken_population - (N)) < 1.0e-08);

    // Repeat this for the RDM of a DOCI expansion


    // Solve the SCF equations
    GQCP::PlainRHFSCFSolver plain_scf_solver (ham_par, N2);  // The DIIS SCF solver seems to find a wrong minimum, so use a plain solver instead
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    ham_par.transform(rhf.get_C());

    GQCP::FockSpace fock_space (K, N/2);
    GQCP::DOCI doci (fock_space);

    GQCP::CISolver ci_solver (doci, ham_par);

    numopt::eigenproblem::DenseSolverOptions solver_options;
    ci_solver.solve(solver_options);

    numopt::eigenproblem::Eigenpair eigen_pair = ci_solver.get_eigenpair(0);

    GQCP::RDMCalculator rdm_calculator (fock_space);

    GQCP::OneRDMs one_rdms = rdm_calculator.calculate1RDMs(eigen_pair.get_eigenvector());

    double mulliken_population_2 = GQCP::calculateExpectationValue(mulliken, one_rdms.one_rdm);
    BOOST_CHECK(std::abs(mulliken_population_2 - (N)) < 1.0e-08);
}
