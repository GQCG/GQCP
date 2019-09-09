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
#define BOOST_TEST_MODULE "HamiltonianParameters"

#include <boost/test/unit_test.hpp>

#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "RHF/PlainRHFSCFSolver.hpp"
#include "Utilities/miscellaneous.hpp"
#include "Utilities/linalg.hpp"

#include <boost/math/constants/constants.hpp>


/*
 *  HELPER FUNCTIONS
 */

/**
 *  @return a toy 2-DM:
 *      d(p,q,r,s) = delta_pq delta_rs
 */
GQCP::TwoRDM<double> calculateToy2DM() {

    GQCP::TwoRDM<double> d (2);
    d.setZero();

    for (size_t p = 0; p < 2; p++) {
        for (size_t q = 0; q < 2; q++) {
            for (size_t r = 0; r < 2; r++) {
                for (size_t s = 0; s < 2; s++) {
                    if ((p == q) && (r == s)) {
                        d(p,q,r,s) = 1;
                    }
                }
            }
        }
    }

    return d;
};



/**
 *  @return toy 2-electron integrals:
 *      g(p,q,r,s) = -0.5 delta_pq delta_rs
 */
GQCP::ScalarSQTwoElectronOperator<double> calculateToyTwoElectronIntegrals() {

    GQCP::QCRankFourTensor<double> g (2);
    g.setZero();

    for (size_t p = 0; p < 2; p++) {
        for (size_t q = 0; q < 2; q++) {
            for (size_t r = 0; r < 2; r++) {
                for (size_t s = 0; s < 2; s++) {
                    if ((p == q) && (r == s)) {
                        g(p,q,r,s) = -0.5;
                    }
                }
            }
        }
    }

    return GQCP::ScalarSQTwoElectronOperator<double>({g});
};



/*
 *  UNIT TESTS - CONSTRUCTORS
 */

BOOST_AUTO_TEST_CASE ( HamiltonianParameters_constructor ) {

    // Create a scalar basis
    auto water = GQCP::Molecule::ReadXYZ("data/h2o.xyz");
    auto ao_basis_ptr = std::make_shared<GQCP::ScalarBasis<GQCP::GTOShell>>(water, "STO-3G");


    // Create One- and SQTwoElectronOperators (and a transformation matrix) with compatible dimensions
    size_t K = ao_basis_ptr->numberOfBasisFunctions();
    GQCP::QCMatrix<double> S = GQCP::QCMatrix<double>::Random(K, K);
    GQCP::QCMatrix<double> H_core = GQCP::QCMatrix<double>::Random(K, K);

    GQCP::QCRankFourTensor<double> g (K);
    g.setRandom();

    GQCP::TransformationMatrix<double> C = GQCP::TransformationMatrix<double>::Random(K, K);


    // Check if a correct constructor works
    GQCP::SQHamiltonian<double> sq_hamiltonian (ao_basis_ptr, GQCP::ScalarSQOneElectronOperator<double>({S}), GQCP::ScalarSQOneElectronOperator<double>({H_core}), GQCP::ScalarSQTwoElectronOperator<double>({g}), C);


    // Check if wrong arguments result in a throw
    GQCP::QCMatrix<double> S_faulty = GQCP::QCMatrix<double>::Random(K+1, K+1);
    GQCP::QCMatrix<double> H_core_faulty = GQCP::QCMatrix<double>::Random(K+1, K+1);

    GQCP::QCRankFourTensor<double> g_faulty (K+1);

    GQCP::TransformationMatrix<double> C_faulty = GQCP::TransformationMatrix<double>::Random(K+1, K+1);

    BOOST_CHECK_THROW(GQCP::SQHamiltonian<double> (ao_basis_ptr, GQCP::ScalarSQOneElectronOperator<double>({S_faulty}), GQCP::ScalarSQOneElectronOperator<double>({H_core}), GQCP::ScalarSQTwoElectronOperator<double>({g}), C), std::invalid_argument);
    BOOST_CHECK_THROW(GQCP::SQHamiltonian<double> (ao_basis_ptr, GQCP::ScalarSQOneElectronOperator<double>({S}), GQCP::ScalarSQOneElectronOperator<double>({H_core_faulty}), GQCP::ScalarSQTwoElectronOperator<double>({g}), C), std::invalid_argument);
    BOOST_CHECK_THROW(GQCP::SQHamiltonian<double> (ao_basis_ptr, GQCP::ScalarSQOneElectronOperator<double>({S}), GQCP::ScalarSQOneElectronOperator<double>({H_core}), GQCP::ScalarSQTwoElectronOperator<double>({g_faulty}), C), std::invalid_argument);
    BOOST_CHECK_THROW(GQCP::SQHamiltonian<double> (ao_basis_ptr, GQCP::ScalarSQOneElectronOperator<double>({S}), GQCP::ScalarSQOneElectronOperator<double>({H_core}), GQCP::ScalarSQTwoElectronOperator<double>({g}), C_faulty), std::invalid_argument);

    // Check if we can't use a zero matrix as overlap matrix
    GQCP::QCMatrix<double> S_zero = GQCP::QCMatrix<double>::Zero(K, K);
    BOOST_CHECK_THROW(GQCP::SQHamiltonian<double> (ao_basis_ptr, GQCP::ScalarSQOneElectronOperator<double>({S_zero}), GQCP::ScalarSQOneElectronOperator<double>({H_core}), GQCP::ScalarSQTwoElectronOperator<double>({g}), C), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( rotate_argument ) {

    // Create a well-behaved Hamiltonian
    size_t K = 3;
    GQCP::QCMatrix<double> S_op = GQCP::QCMatrix<double>::Random(K, K);
    GQCP::QCMatrix<double> H_op = GQCP::QCMatrix<double>::Random(K, K);
    GQCP::QCRankFourTensor<double> g_op (K);
    g_op.setRandom();

    GQCP::SQHamiltonian<double> sq_hamiltonian (nullptr, GQCP::ScalarSQOneElectronOperator<double>({S_op}), GQCP::ScalarSQOneElectronOperator<double>({H_op}), GQCP::ScalarSQTwoElectronOperator<double>({g_op}), GQCP::TransformationMatrix<double>::Random(K, K));


    // Check if we can't rotate with a non-unitary matrix
    GQCP::TransformationMatrix<double> T (K);
    T << 0.5, 0.5, -2.0,
         3.0, 0.0,  1.5,
         0.0, 0.0,  2.5;
    BOOST_CHECK_THROW(sq_hamiltonian.rotate(GQCP::TransformationMatrix<double>(T)), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( rotate_overlap_matrix ) {

    // Check if a rotation that interchanges two orbitals changes the non-identity overlap matrix
    GQCP::JacobiRotationParameters jacobi_rotation_parameters {1, 0, boost::math::constants::half_pi<double>()};  // interchanges two orbitals

    size_t K = 3;
    GQCP::QCMatrix<double> S_op (K);
    S_op << 1.0, 0.5, 0.0,
            0.5, 2.0, 0.0,
            0.0, 0.0, 1.0;

    GQCP::QCMatrix<double> S_rotated_ref (K);  // manual calculation
    S_rotated_ref <<  2.0, -0.5, 0.0,
                     -0.5,  1.0, 0.0,
                      0.0,  0.0, 1.0;

    GQCP::QCMatrix<double> H_op = GQCP::QCMatrix<double>::Random(K, K);

    GQCP::QCRankFourTensor<double> g_op (K);
    g_op.setRandom();


    // Check the Jacobi rotation
    GQCP::SQHamiltonian<double> ham_par_jacobi (nullptr, GQCP::ScalarSQOneElectronOperator<double>({S_op}), GQCP::ScalarSQOneElectronOperator<double>({H_op}), GQCP::ScalarSQTwoElectronOperator<double>({g_op}), GQCP::TransformationMatrix<double>::Random(K, K));
    ham_par_jacobi.rotate(jacobi_rotation_parameters);
    BOOST_CHECK(ham_par_jacobi.get_S().parameters().isApprox(S_rotated_ref, 1.0e-08));


    // Check for a unitary transformation
    GQCP::SQHamiltonian<double> sq_hamiltonian (nullptr, GQCP::ScalarSQOneElectronOperator<double>({S_op}), GQCP::ScalarSQOneElectronOperator<double>({H_op}), GQCP::ScalarSQTwoElectronOperator<double>({g_op}), GQCP::TransformationMatrix<double>::Random(K, K));
    auto J = GQCP::TransformationMatrix<double>::FromJacobi(jacobi_rotation_parameters, K);
    sq_hamiltonian.rotate(GQCP::TransformationMatrix<double>(J));
    BOOST_CHECK(sq_hamiltonian.get_S().parameters().isApprox(S_rotated_ref, 1.0e-08));
}


BOOST_AUTO_TEST_CASE ( constructor_C ) {

    // Create a dummy Hamiltonian
    std::shared_ptr<GQCP::ScalarBasis<GQCP::GTOShell>> ao_basis;
    size_t K = 4;
    GQCP::QCMatrix<double> S = GQCP::QCMatrix<double>::Random(K, K);
    GQCP::QCMatrix<double> H_core = GQCP::QCMatrix<double>::Random(K, K);

    GQCP::QCRankFourTensor<double> g (K);

    GQCP::TransformationMatrix<double> C = GQCP::TransformationMatrix<double>::Random(K, K);

    GQCP::SQHamiltonian<double> sq_hamiltonian (ao_basis, GQCP::ScalarSQOneElectronOperator<double>({S}), GQCP::ScalarSQOneElectronOperator<double>({H_core}), GQCP::ScalarSQTwoElectronOperator<double>({g}), C);


    // Check if we can create a transformed Hamiltonian
    GQCP::TransformationMatrix<double> T = GQCP::TransformationMatrix<double>::Random(K, K);
    GQCP::SQHamiltonian<double> transformed_random_hamiltonian_parameters (sq_hamiltonian, T);
}


/*
 *  UNIT TESTS - NAMED CONSTRUCTORS
 */

BOOST_AUTO_TEST_CASE ( constructMolecularHamiltonianParameters ) {

    // Set up a basis
    auto h2 = GQCP::Molecule::ReadXYZ("data/h2_szabo.xyz");


    // Check if we can construct the molecular Hamiltonian
    GQCP::SingleParticleBasis<double, GQCP::GTOShell> sp_basis (h2, "STO-3G");
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(sp_basis, h2);  // in an AO basis
    auto g = sq_hamiltonian.twoElectron().parameters();


    // Check with reference values from Szabo
    GQCP::QCMatrix<double> ref_S (2);
    ref_S << 1.0,    0.6593,
             0.6593, 1.0;

    GQCP::QCMatrix<double> ref_H_core (2);
    ref_H_core << -1.1204, -0.9584,
                  -0.9584, -1.1204;


    BOOST_CHECK(sq_hamiltonian.get_S().parameters().isApprox(ref_S, 1.0e-04));
    BOOST_CHECK(sq_hamiltonian.core().parameters().isApprox(ref_H_core, 1.0e-04));

    BOOST_CHECK(std::abs(g(0,0,0,0) - 0.7746) < 1.0e-04);
    BOOST_CHECK(std::abs(g(0,0,0,0) - g(1,1,1,1)) < 1.0e-12);
    BOOST_CHECK(std::abs(g(0,0,1,1) - 0.5697) < 1.0e-04);
    BOOST_CHECK(std::abs(g(1,0,0,0) - 0.4441) < 1.0e-04);
    BOOST_CHECK(std::abs(g(1,0,0,0) - g(1,1,1,0)) < 1.0e-12);
    BOOST_CHECK(std::abs(g(1,0,1,0) - 0.2970) < 1.0e-04);
}


BOOST_AUTO_TEST_CASE ( FCIDUMP_reader ) {

    auto fcidump_ham_par = GQCP::SQHamiltonian<double>::ReadFCIDUMP("data/beh_cation_631g_caitlin.FCIDUMP");

    // Check if the one-electron integrals are read in correctly from a previous implementation
    GQCP::QCMatrix<double> h_SO = fcidump_ham_par.core().parameters();

    BOOST_CHECK(std::abs(h_SO(0,0) - (-8.34082)) < 1.0e-5);
    BOOST_CHECK(std::abs(h_SO(5,1) - 0.381418) < 1.0e-6);
    BOOST_CHECK(std::abs(h_SO(14,0) - 0.163205) < 1.0e-6);
    BOOST_CHECK(std::abs(h_SO(13,6) - (-5.53204e-16)) < 1.0e-16);
    BOOST_CHECK(std::abs(h_SO(15,11) - (-0.110721)) < 1.0e-6);


    // Check if the two-electron integrals are read in correctly from a previous implementation
    GQCP::QCRankFourTensor<double> g_SO = fcidump_ham_par.twoElectron().parameters();

    BOOST_CHECK(std::abs(g_SO(2,5,4,4) - 0.0139645) < 1.0e-6);
    BOOST_CHECK(std::abs(g_SO(2,6,3,0) - 5.16622e-18) < 1.0e-17);
    BOOST_CHECK(std::abs(g_SO(3,1,3,0) - (-0.0141251)) <  1.0e-6);
    BOOST_CHECK(std::abs(g_SO(4,6,4,6) - 0.0107791) < 1.0e-6);
    BOOST_CHECK(std::abs(g_SO(4,15,11,1) - (9.33375e-19)) < 1.0e-17);
    BOOST_CHECK(std::abs(g_SO(6,10,5,9) - (-3.81422e-18)) < 1.0e-17);
    BOOST_CHECK(std::abs(g_SO(7,7,2,1) - (-0.031278)) < 1.0e-6);
    BOOST_CHECK(std::abs(g_SO(8,15,9,9) - (-2.80093e-17)) < 1.0e-16);
    BOOST_CHECK(std::abs(g_SO(9,14,0,9) - 0.00161985) < 1.0e-7);
    BOOST_CHECK(std::abs(g_SO(10,1,4,3) - 0.00264603) < 1.0e-7);
    BOOST_CHECK(std::abs(g_SO(11,4,9,3) - (-0.0256623)) < 1.0e-6);
    BOOST_CHECK(std::abs(g_SO(12,9,0,4) - 0.0055472) < 1.0e-6);
    BOOST_CHECK(std::abs(g_SO(13,15,15,13) - 0.00766898) < 1.0e-7);
    BOOST_CHECK(std::abs(g_SO(14,2,12,3) - 0.0104266) < 1.0e-7);
    BOOST_CHECK(std::abs(g_SO(15,5,10,10) - 0.00562608) < 1.0e-7);
}


BOOST_AUTO_TEST_CASE ( FCIDUMP_reader_HORTON ) {

    // Check the same reference value that HORTON does
    auto fcidump_ham_par = GQCP::SQHamiltonian<double>::ReadFCIDUMP("data/h2_psi4_horton.FCIDUMP");

    GQCP::QCRankFourTensor<double> g_SO = fcidump_ham_par.twoElectron().parameters();
    BOOST_CHECK(std::abs(g_SO(6,5,1,0) - 0.0533584656) <  1.0e-7);
}



/*
 *  UNIT TESTS - METHODS
 */

BOOST_AUTO_TEST_CASE ( calculate_generalized_Fock_matrix_and_super_invalid_arguments ) {

    // Initialize toy HamiltonianParameters
    std::shared_ptr<GQCP::ScalarBasis<GQCP::GTOShell>> ao_basis;
    GQCP::QCMatrix<double> S = GQCP::QCMatrix<double>::Identity(2, 2);
    GQCP::QCMatrix<double> h = GQCP::QCMatrix<double>::Zero(2, 2);
    GQCP::QCRankFourTensor<double> g (2);
    GQCP::SQHamiltonian<double> sq_hamiltonian (ao_basis, GQCP::ScalarSQOneElectronOperator<double>({S}), GQCP::ScalarSQOneElectronOperator<double>({h}), GQCP::ScalarSQTwoElectronOperator<double>({g}), GQCP::TransformationMatrix<double>::Identity(2, 2));


    // Create valid and invalid density matrices (with respect to the dimensions of the SOBasis)
    GQCP::OneRDM<double> D_valid = GQCP::OneRDM<double>::Zero(2, 2);
    GQCP::OneRDM<double> D_invalid = GQCP::OneRDM<double>::Zero(3, 3);

    GQCP::TwoRDM<double> d_valid (2);
    GQCP::TwoRDM<double> d_invalid (3);



    // Test a faulty function calls
    BOOST_REQUIRE_THROW(sq_hamiltonian.calculateFockianMatrix(D_invalid, d_valid), std::invalid_argument);
    BOOST_REQUIRE_THROW(sq_hamiltonian.calculateFockianMatrix(D_valid, d_invalid), std::invalid_argument);

    BOOST_REQUIRE_THROW(sq_hamiltonian.calculateSuperFockianMatrix(D_invalid, d_valid), std::invalid_argument);
    BOOST_REQUIRE_THROW(sq_hamiltonian.calculateSuperFockianMatrix(D_valid, d_invalid), std::invalid_argument);


    // Test correct function calls
    sq_hamiltonian.calculateFockianMatrix(D_valid, d_valid);
    sq_hamiltonian.calculateSuperFockianMatrix(D_valid, d_valid);
}


BOOST_AUTO_TEST_CASE ( calculate_Fockian_and_super ) {

    // We test the function by a manual calculation of nonsensical toy 1- and 2-RDMS and one- and two-electron integrals
    // Set up the toy 1- and 2-RDMs
    GQCP::OneRDM<double> D (2);
    D << 0, 1,
         1, 0;

    auto d = calculateToy2DM();

    // Set up the toy Hamiltonian
    std::shared_ptr<GQCP::ScalarBasis<GQCP::GTOShell>> ao_basis;
    GQCP::QCMatrix<double> S = GQCP::QCMatrix<double>::Identity(2, 2);
    GQCP::QCMatrix<double> h = GQCP::QCMatrix<double>::Zero(2, 2);
    h << 1, 0,
         0, 1;

    auto g = calculateToyTwoElectronIntegrals();
    GQCP::SQHamiltonian<double> sq_hamiltonian (ao_basis, GQCP::ScalarSQOneElectronOperator<double>({S}), GQCP::ScalarSQOneElectronOperator<double>({h}), g, GQCP::TransformationMatrix<double>::Identity(2, 2));


    // Construct the reference Fockian matrix
    GQCP::SquareMatrix<double> F_ref (2);
    F_ref << -1.00,  1.00,
              1.00, -1.00;

    // Construct the reference super generalized Fock matrix
    GQCP::QCRankFourTensor<double> G_ref (2);
    G_ref.setZero();
    for (size_t p = 0; p < 2; p++) {
        for (size_t q = 0; q < 2; q++) {
            for (size_t r = 0; r < 2; r++) {
                for (size_t s = 0; s < 2; s++) {
                    if (q == r) {
                        G_ref(p,q,r,s) += F_ref(p,s);
                    }

                    // One-electron part is simplified by manual calculation
                    if (p == s) {
                        G_ref(p,q,r,s) -= D(r,q);
                    }

                    // Two-electron part is simplified by manual calculation
                    if ((p == s) && (q == r)) {
                        G_ref(p,q,r,s) += 1;
                    }
                }
            }
        }
    }


    BOOST_CHECK(F_ref.isApprox(sq_hamiltonian.calculateFockianMatrix(D, d), 1.0e-12));
    BOOST_CHECK(G_ref.isApprox(sq_hamiltonian.calculateSuperFockianMatrix(D, d), 1.0e-12));
}


BOOST_AUTO_TEST_CASE ( calculateEdmistonRuedenbergLocalizationIndex ) {

    // Create a toy Hamiltonian: only the two-electron integrals are important
    size_t K = 5;
    GQCP::QCMatrix<double> S_op = GQCP::QCMatrix<double>::Identity(K, K);
    GQCP::QCMatrix<double> H_op = GQCP::QCMatrix<double>::Random(K, K);

    GQCP::QCRankFourTensor<double> g_op (K);
    g_op.setZero();
    for (size_t i = 0; i < K; i++) {
        g_op(i,i,i,i) = 2*static_cast<float>(i);
    }

    GQCP::SQHamiltonian<double> sq_hamiltonian (nullptr, GQCP::ScalarSQOneElectronOperator<double>({S_op}), GQCP::ScalarSQOneElectronOperator<double>({H_op}), GQCP::ScalarSQTwoElectronOperator<double>({g_op}), GQCP::TransformationMatrix<double>::Identity(K, K));


    BOOST_CHECK(std::abs(sq_hamiltonian.calculateEdmistonRuedenbergLocalizationIndex(3) - 6.0) < 1.0e-08);
    BOOST_CHECK(std::abs(sq_hamiltonian.calculateEdmistonRuedenbergLocalizationIndex(4) - 12.0) < 1.0e-08);
}


BOOST_AUTO_TEST_CASE ( areOrbitalsOrthonormal ) {

    // We assume that the orbitals in an FCIDUMP file are orthonormal
    auto ham_par_fcidump = GQCP::SQHamiltonian<double>::ReadFCIDUMP("data/h2_psi4_horton.FCIDUMP");
    BOOST_CHECK(ham_par_fcidump.areOrbitalsOrthonormal());


    // The orbitals in an AO basis are not orthonormal
    auto h2o = GQCP::Molecule::ReadXYZ("data/h2o.xyz");
    GQCP::SingleParticleBasis<double, GQCP::GTOShell> sp_basis (h2o, "STO-3G");
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(sp_basis, h2o);  // in an AO basis
    BOOST_CHECK(!sq_hamiltonian.areOrbitalsOrthonormal());


    // The orbitals in the RHF basis are orthonormal
    GQCP::PlainRHFSCFSolver plain_scf_solver (sq_hamiltonian, h2o);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();
    sq_hamiltonian.transform(rhf.get_C());
    BOOST_CHECK(sq_hamiltonian.areOrbitalsOrthonormal());
}


BOOST_AUTO_TEST_CASE ( dissociatedMoleculeParameters ) {

    // Test if we can succesfully initialize NO+ at long intra molecular distance
    auto N = GQCP::Nucleus(7, 3.5, 0, 0);
    auto O = GQCP::Nucleus(8, -3.5, 0, 0);
    std::vector<GQCP::Nucleus> nuclei {N,O};
    auto NO = GQCP::Molecule(nuclei, +1);

    GQCP::SingleParticleBasis<double, GQCP::GTOShell> sp_basis (NO, "STO-3G");
    BOOST_CHECK_NO_THROW(GQCP::SQHamiltonian<double>::Molecular(sp_basis, NO));
}
