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
#define BOOST_TEST_MODULE "expectation_values"

#include <boost/test/unit_test.hpp>

#include "Properties/expectation_values.hpp"

#include "Basis/Integrals/Interfaces/LibintInterfacer.hpp"
#include "CISolver/CISolver.hpp"
#include "HamiltonianBuilder/DOCI.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "RDM/RDMCalculator.hpp"
#include "RHF/DIISRHFSCFSolver.hpp"
#include "RHF/PlainRHFSCFSolver.hpp"
#include "units.hpp"


BOOST_AUTO_TEST_CASE ( one_electron_throw ) {

    GQCP::ScalarSQOneElectronOperator<double> h ({GQCP::QCMatrix<double>::Zero(2, 2)});
    GQCP::OneRDM<double> D_valid = GQCP::OneRDM<double>::Zero(2, 2);
    GQCP::OneRDM<double> D_invalid = GQCP::OneRDM<double>::Zero(3, 3);

    BOOST_CHECK_THROW(GQCP::calculateExpectationValue(h, D_invalid), std::invalid_argument);
    BOOST_CHECK_NO_THROW(GQCP::calculateExpectationValue(h, D_valid));
}


BOOST_AUTO_TEST_CASE ( two_electron_throw ) {

    GQCP::ScalarSQTwoElectronOperator<double> g (2);

    GQCP::TwoRDM<double> d_valid (2);
    GQCP::TwoRDM<double> d_invalid (3);

    BOOST_CHECK_THROW(GQCP::calculateExpectationValue(g, d_invalid), std::invalid_argument);
    BOOST_CHECK_NO_THROW(GQCP::calculateExpectationValue(g, d_valid));
}


BOOST_AUTO_TEST_CASE ( mulliken_N2_STO_3G ) {

    // Check that the mulliken population of N2 is 14 (N)

    // Initialize the molecule and the molecular Hamiltonian for N2
    GQCP::Nucleus N_1 (7, 0.0, 0.0, 0.0);
    GQCP::Nucleus N_2 (7, 0.0, 0.0, GQCP::units::angstrom_to_bohr(1.134));  // from CCCBDB, STO-3G geometry
    std::vector<GQCP::Nucleus> nuclei {N_1, N_2};
    GQCP::Molecule N2 (nuclei);

    GQCP::SingleParticleBasis<double, GQCP::GTOShell> sp_basis (N2, "STO-3G");
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(sp_basis, N2);  // in an AO basis
    size_t K = sq_hamiltonian.get_K();

    // We include all basis functions
    GQCP::Vectoru gto_list (K);
    for(size_t i = 0; i<K; i++){
        gto_list[i] = i;
    }

    GQCP::ScalarSQOneElectronOperator<double> mulliken = sp_basis.calculateMullikenOperator(gto_list);

    size_t N = N2.numberOfElectrons();

    // Create a 1-RDM for N2
    GQCP::OneRDM<double> one_rdm = GQCP::calculateRHF1RDM(K, N);

    double mulliken_population = GQCP::calculateExpectationValue(mulliken, one_rdm)[0];
    BOOST_CHECK(std::abs(mulliken_population - (N)) < 1.0e-06);


    // Repeat this for a DOCI-RDM

    // Solve the SCF equations
    GQCP::PlainRHFSCFSolver plain_scf_solver (sq_hamiltonian, N2);  // the DIIS SCF solver seems to find a wrong minimum, so use a plain solver instead
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    sq_hamiltonian.transform(rhf.get_C());

    GQCP::FockSpace fock_space (K, N/2);
    GQCP::DOCI doci (fock_space);

    GQCP::CISolver ci_solver (doci, sq_hamiltonian);

    GQCP::DenseSolverOptions solver_options;
    ci_solver.solve(solver_options);

    GQCP::RDMCalculator rdm_calculator (ci_solver.makeWavefunction());

    GQCP::OneRDMs<double> one_rdms = rdm_calculator.calculate1RDMs();

    double mulliken_population_2 = GQCP::calculateExpectationValue(mulliken, one_rdms.one_rdm)[0];
    BOOST_CHECK(std::abs(mulliken_population_2 - (N)) < 1.0e-06);
}
