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
#define BOOST_TEST_MODULE "constrained_RHF"


#include "RHF/DIISRHFSCFSolver.hpp"
#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "properties/expectation_values.hpp"
#include <random>

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


BOOST_AUTO_TEST_CASE ( constrained_CO_test ) {

    // Create a Molecule and an AOBasis with the assumed geometry
    auto CO = GQCP::Molecule::Readxyz("../tests/data/CO_mulliken.xyz");
    auto ao_ham_par = GQCP::HamiltonianParameters::Molecular(CO, "STO-3G");

    size_t K = ao_ham_par.get_K();
    size_t N = CO.get_N();

    GQCP::OneRDM one_rdm = GQCP::calculateRHF1RDM(K, N);

    // Initialize the ref data form:
    // "Self-consistent methods constrained to a fixed number of particles in a given fragment
    // and its relation to the electronegativity equalization method"
    // Andrés Cedillo • Dimitri Van Neck • Patrick Bultinck
    // DOI 10.1007/s00214-012-1227-6
    Eigen::Matrix<double, 21, 3> CO_data;
    //          lambda  charge  energy
    CO_data <<  -1.0,  1.73, -110.530475,
                -0.9,  1.62, -110.634766,
                -0.8,  1.50, -110.737606,
                -0.7,  1.37, -110.836596,
                -0.6,  1.23, -110.929219,
                -0.5,  1.08, -111.012983,
                -0.4,  0.91, -111.085538,
                -0.3,  0.75, -111.144760,
                -0.2,  0.57, -111.188802,
                -0.1,  0.39, -111.216114,
                 0,    0.20, -111.225446,
                 0.1,  0.01, -111.215842,
                 0.2, -0.19, -111.186619,
                 0.3, -0.38, -111.137356,
                 0.4, -0.58, -111.067872,
                 0.5, -0.78, -110.978210,
                 0.6, -0.98, -110.868621,
                 0.7, -1.18, -110.739544,
                 0.8, -1.38, -110.591599,
                 0.9, -1.57, -110.425574,
                 1.0, -1.77, -110.242423;

    // Pick a set of AO's (GTOs of the carbon atom)
    GQCP::Vectoru gto_list = {0,1,2,3,4};

    // Iterate over multipliers
    size_t index = 0;
    for (double i = -1.0; i < 1.01; i += 0.1) {
        // Calculate the Mulliken operator
        auto mulliken_operator = ao_ham_par.calculateMullikenOperator(gto_list);

        // Contrain the original Hamiltonian parameters
        auto constrained_ham_par = ao_ham_par.constrain(mulliken_operator, i);

        // Create a DIIS RHF SCF solver and solve the SCF equations
        GQCP::DIISRHFSCFSolver diis_scf_solver (constrained_ham_par, CO);
        diis_scf_solver.solve();
        auto rhf = diis_scf_solver.get_solution();

        // Transform the mulliken operator to the basis in which the RHF energies are calculated
        mulliken_operator.transform(rhf.get_C());

        // Retrieve the RHF "energy"
        double expectation_value = rhf.get_electronic_energy();

        // Retrieve the expectation value of the Mulliken operator (aka the population)
        double mulliken_population = GQCP::calculateExpectationValue(mulliken_operator, one_rdm);

        // Retrieve the total energy by adding the lambda times the expectation value of the constraining operator
        double total_energy = expectation_value + i * mulliken_population + CO.calculateInternuclearRepulsionEnergy();

        // Mulliken charge on the carbon atom
        double C_charge = 6 - mulliken_population;

        // The lenient threshold is because of we are not certain of exact the geometry in the paper
        BOOST_CHECK(std::abs(total_energy - CO_data(index, 2)) < 1.0e-4);
        BOOST_CHECK(std::abs(C_charge - CO_data(index, 1)) < 1.0e-2);

        index++;
    }
}


BOOST_AUTO_TEST_CASE ( constrained_CO_test_random_transformation) {
    // Repeat the same test but perform a random transformation
    // The Hartree-Fock procedure should be invariant under random transformations
    // This tests if our Mulliken operator remains correct if we transform before the procedure.

    // Create a Molecule and an AOBasis with the assumed geometry
    auto CO = GQCP::Molecule::Readxyz("../tests/data/CO_mulliken.xyz");
    auto ao_ham_par = GQCP::HamiltonianParameters::Molecular(CO, "STO-3G");

    size_t K = ao_ham_par.get_K();
    size_t N = CO.get_N();

    Eigen::MatrixXd T = Eigen::MatrixXd::Random(K, K);
    // set diagonal elements to 1
    for (int i = 0; i < K; i++) {
        T(i,i) = 1;
    }

    ao_ham_par.transform(T);

    GQCP::OneRDM one_rdm = GQCP::calculateRHF1RDM(K, N);

    // Initialize the ref data form:
    // Self-consistent methods constrained to a fixed number of particles in a given fragment
    // and its relation to the electronegativity equalization method
    // Andrés Cedillo • Dimitri Van Neck • Patrick Bultinck
    // DOI 10.1007/s00214-012-1227-6
    Eigen::Matrix<double,21,3> CO_data;
    //          lambda  charge  energy
    CO_data <<  -1.0 , 1.73  , -110.530475 ,
                -0.9 , 1.62  , -110.634766 ,
                -0.8 , 1.50  , -110.737606 ,
                -0.7 , 1.37  , -110.836596 ,
                -0.6 , 1.23  , -110.929219 ,
                -0.5 , 1.08  , -111.012983 ,
                -0.4 , 0.91  , -111.085538 ,
                -0.3 , 0.75  , -111.144760 ,
                -0.2 , 0.57  , -111.188802 ,
                -0.1 , 0.39  , -111.216114 ,
                0    , 0.20  , -111.225446 ,
                0.1  , 0.01  , -111.215842 ,
                0.2  , -0.19 , -111.186619 ,
                0.3  , -0.38 , -111.137356 ,
                0.4  , -0.58 , -111.067872 ,
                0.5  , -0.78 , -110.978210 ,
                0.6  , -0.98 , -110.868621 ,
                0.7  , -1.18 , -110.739544 ,
                0.8  , -1.38 , -110.591599 ,
                0.9  , -1.57 , -110.425574 ,
                1.0  , -1.77 , -110.242423 ;

    // Pick a set of AO's (GTOs of the carbon atom)
    GQCP::Vectoru gto_list = {0,1,2,3,4};

    // Iterate over multipliers
    size_t index = 0;
    for (double i = -1; i<1.01; i += 0.1) {
        // Calculate the Mulliken operator
        auto mulliken_operator = ao_ham_par.calculateMullikenOperator(gto_list);

        // Contrain the original Hamiltonian parameters
        auto constrained_ham_par = ao_ham_par.constrain(mulliken_operator, i);

        // Create a DIIS RHF SCF solver and solve the SCF equations
        GQCP::DIISRHFSCFSolver diis_scf_solver (constrained_ham_par, CO);
        diis_scf_solver.solve();
        auto rhf = diis_scf_solver.get_solution();

        // Transform the mulliken operator to the basis in which the RHF energies are calculated
        mulliken_operator.transform(rhf.get_C());

        // Retrieve the RHF "energy"
        double expectation_value = rhf.get_electronic_energy();

        // Retrieve the expectation value of the Mulliken operator (aka the population)
        double mulliken_population = GQCP::calculateExpectationValue(mulliken_operator, one_rdm);

        // Retrieve the total energy by adding the lambda times the expectation value of the constraining operator
        double total_energy = expectation_value + i * mulliken_population + CO.calculateInternuclearRepulsionEnergy();

        // Mulliken charge on the carbon atom
        double C_charge = 6 - mulliken_population;

        // The lenient threshold is because of we are not certain of exact the geometry in the paper
        BOOST_CHECK(std::abs(total_energy - CO_data(index, 2)) < 1.0e-4);
        BOOST_CHECK(std::abs(C_charge - CO_data(index, 1)) < 1.0e-2);

        index++;
    }
}
