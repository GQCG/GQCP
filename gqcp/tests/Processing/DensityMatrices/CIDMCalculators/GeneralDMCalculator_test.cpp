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

#define BOOST_TEST_MODULE "RDMCalculator_test"

#include <boost/test/unit_test.hpp>

#include "ONVBasis/SpinResolvedONVBasis.hpp"
#include "ONVBasis/SpinUnresolvedONVBasis.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "Processing/DensityMatrices/CIDMCalculators/GeneralDMCalculator.hpp"
#include "Processing/DensityMatrices/CIDMCalculators/GeneralSpinUnresolvedDMCalculator.hpp"
#include "QCMethod/CI/HamiltonianBuilder/DOCI.hpp"


// BOOST_AUTO_TEST_CASE(constructor) {

//     // Test polymorphic entry for RDM (from SeniorityZeroDMCalculator test-case).

//     // Get the 1-RDM from DOCI
//     size_t N = 4;  // 4 electrons
//     auto sq_hamiltonian = GQCP::SQHamiltonian<double>::ReadFCIDUMP("data/lih_631g_caitlin.FCIDUMP");
//     size_t K = sq_hamiltonian.dimension();  // 16 SO

//     // Abstract pointer to test RDM
//     std::shared_ptr<GQCP::BaseONVBasis> fock_space_dy {new GQCP::SpinUnresolvedONVBasis(K, N / 2)};  // dim = 120
//     GQCP::SpinUnresolvedONVBasis fock_space {K, N / 2};                                              // dim = 120

//     GQCP::DOCI doci {fock_space};

//     // Specify solver options and solve the eigenvalue problem
//     // Solve the dense DOCI eigenvalue problem
//     GQCP::CISolver ci_solver {doci, sq_hamiltonian};
//     GQCP::DenseSolverOptions solver_options;
//     ci_solver.solve(solver_options);

//     GQCP::VectorX<double> coef = ci_solver.eigenpair().eigenvector();

//     // Check if the DOCI 1-RDM has the proper trace.
//     GQCP::GeneralDMCalculator doci_rdm {*fock_space_dy};
//     doci_rdm.setCoefficients(coef);
//     GQCP::SpinResolvedOneDM<double> one_rdms = doci_rdm.calculate1RDMs();

//     BOOST_CHECK(std::abs(one_rdms.spinSummed().trace() - N) < 1.0e-12);
// }

BOOST_AUTO_TEST_CASE(no_vector_throws) {

    size_t K = 5;
    size_t N = 4;
    GQCP::SpinUnresolvedONVBasis fock_space {K, N};

    GQCP::VectorX<double> coeff {fock_space.dimension()};
    coeff << 1, 1, -2, 4, -5;

    // Test if throws when no vector is set
    GQCP::GeneralDMCalculator doci_rdm {fock_space};
    BOOST_CHECK_THROW(doci_rdm.calculate1RDMs(), std::logic_error);
}


BOOST_AUTO_TEST_CASE(operator_call_throw) {

    // Create a test wave function
    size_t M = 3;
    size_t N = 1;
    GQCP::SpinUnresolvedONVBasis fock_space {M, N};

    GQCP::VectorX<double> coeff(fock_space.dimension());
    coeff << 1, 2, -3;
    GQCP::GeneralSpinUnresolvedDMCalculator d {fock_space};
    d.setCoefficients(coeff);
    BOOST_CHECK_THROW(d(0), std::invalid_argument);  // need an even number of indices
}


BOOST_AUTO_TEST_CASE(operator_call) {

    // Create a test wave function
    size_t M = 3;
    size_t N = 2;
    GQCP::SpinUnresolvedONVBasis fock_space {M, N};

    GQCP::VectorX<double> coeff(fock_space.dimension());
    coeff << 1, 2, -3;
    GQCP::GeneralSpinUnresolvedDMCalculator d {fock_space};
    d.setCoefficients(coeff);

    BOOST_CHECK(std::abs(d(0, 1, 1, 2) - (-3.0)) < 1.0e-12);  // d(0,1,1,2) : a^\dagger_0 a^\dagger_1 a_2 a_1
    BOOST_CHECK(std::abs(d(2, 0, 0, 1) - (-2.0)) < 1.0e-12);  // d(2,0,0,1) : a^\dagger_2 a^\dagger_0 a^1 a_0
    BOOST_CHECK(std::abs(d(0, 2, 2, 0) - (-4.0)) < 1.0e-12);  // d(0,2,2,0) : a^\dagger_0 a^dagger_2 a_0 a_2
    BOOST_CHECK(std::abs(d(0, 2, 0, 0) - 0.0) < 1.0e-12);     // d(0,2,0,0) : a^\dagger_0 a^dagger_0 a_0 a_2, double annihilation gives 0.0
}
