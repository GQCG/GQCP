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

#define BOOST_TEST_MODULE "SimpleSQTwoElectronOperator"

#include <boost/test/unit_test.hpp>

#include "Operator/SecondQuantized/RSQTwoElectronOperator.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "Utilities/miscellaneous.hpp"

#include <boost/math/constants/constants.hpp>


/**
 *  In this test suite, we'll test the behavior of `SimpleSQTwoElectronOperator` through a derived class `RSQTwoElectronOperator`.
 * 
 *  We don't have to re-test any methods that are enabled by conforming to `VectorSpaceArithmetic`, since this has been tested in `SimpleSQOneElectronOperator` which also derives from `SQOperatorStorage`.
 */


/*
 *  MARK: Helper functions
 */


/*
 *  Set up toy two-electron integrals, where
 *      g(i,j,k,l) = (i + 1) + 2 * (j + 1) + 4 * (k + 1) + 8 * (l + 1)
 * 
 *  @param dim      The dimension of each axis.
 */
GQCP::ScalarRSQTwoElectronOperator<double> toyTwoElectronIntegrals(const size_t dim) {

    // Initialize a tensor and convert it to an operator.
    auto g = GQCP::SquareRankFourTensor<double>::Zero(dim);

    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
            for (size_t k = 0; k < dim; k++) {
                for (size_t l = 0; l < dim; l++) {
                    g(i, j, k, l) = (i + 1) + 2 * (j + 1) + 4 * (k + 1) + 8 * (l + 1);
                }
            }
        }
    }

    return GQCP::ScalarRSQTwoElectronOperator<double> {g};
}


/*
 *  Set up toy two-electron integrals, where
 *      g(i,j,k,l) = 1
 * 
 *  @param dim      The dimension of each axis.
 */
GQCP::ScalarRSQTwoElectronOperator<double> toyTwoElectronIntegrals2(const size_t dim) {

    // Initialize a tensor and convert it into an operator.
    auto g = GQCP::SquareRankFourTensor<double>::Zero(dim);

    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
            for (size_t k = 0; k < dim; k++) {
                for (size_t l = 0; l < dim; l++) {
                    g(i, j, k, l) = 1;
                }
            }
        }
    }

    return GQCP::ScalarRSQTwoElectronOperator<double> {g};
}


/*
 *  MARK: Tests
 */

/**
 *  Check the interface for constructing SQTwoElectronOperators from `SquareRankFourTensor`s.
 */
BOOST_AUTO_TEST_CASE(constructor) {

    // Check a correct constructor.
    const GQCP::SquareRankFourTensor<double> tensor {3};
    BOOST_CHECK_NO_THROW(GQCP::ScalarRSQTwoElectronOperator<double> O {tensor});


    // Check a faulty constructor.
    const GQCP::Tensor<double, 4> tensor2 {3, 3, 3, 2};
    BOOST_CHECK_THROW(GQCP::ScalarRSQTwoElectronOperator<double> O2 {tensor2}, std::invalid_argument);  // The last axis has an unexpected dimension.
}


/**
 *  Check if the `Zero` named constructor really sets its parameters to all zeros.
 */
BOOST_AUTO_TEST_CASE(Zero) {

    const size_t dim = 2;
    const auto op = GQCP::ScalarRSQTwoElectronOperator<double>::Zero(dim);

    // Create a reference zero tensor.
    GQCP::SquareRankFourTensor<double> ref = GQCP::SquareRankFourTensor<double>::Zero(dim);

    BOOST_CHECK_EQUAL(op.numberOfOrbitals(), dim);
    BOOST_CHECK(op.parameters().isApprox(ref.setZero(), 1.0e-08));
}


/**
 *  Check if the formulas in effectiveOneElectronPartition are implemented correctly.
 */
BOOST_AUTO_TEST_CASE(effectiveOneElectronPartition) {

    const size_t K = 4;
    const auto K_ = static_cast<double>(K);  // Needed in the manual calculation.

    // Set up toy 2-electron integrals.
    const auto g = toyTwoElectronIntegrals(K);

    // Set up the reference effective one-electron integrals by manual calculation.
    auto k_par_ref = GQCP::SquareMatrix<double>::Zero(K);  // 'par_ref' for reference parameters
    for (size_t p = 0; p < K; p++) {
        for (size_t q = 0; q < K; q++) {
            auto p_ = static_cast<double>(p) + 1;
            auto q_ = static_cast<double>(q) + 1;

            k_par_ref(p, q) = -K_ / 2 * (p_ + 8 * q_ + 3 * K_ + 3);
        }
    }

    BOOST_CHECK(k_par_ref.isApprox(g.effectiveOneElectronPartition().parameters(), 1.0e-08));
}


/**
 *  Check if calculateExpectationValue throws when necessary.
 */
BOOST_AUTO_TEST_CASE(calculateExpectationValue_throw) {

    const GQCP::ScalarRSQTwoElectronOperator<double> g {2};

    const GQCP::Orbital2DM<double> d_valid {2};
    const GQCP::Orbital2DM<double> d_invalid {3};

    BOOST_CHECK_THROW(g.calculateExpectationValue(d_invalid), std::invalid_argument);
    BOOST_CHECK_NO_THROW(g.calculateExpectationValue(d_valid));
}


/**
 *  Check whether or not calculateExpectationValue shows the correct behaviour
 */
BOOST_AUTO_TEST_CASE(calculateExpectationValue_behaviour) {

    const size_t dim = 2;

    // Set up toy two-electron integrals.
    const auto op = toyTwoElectronIntegrals(dim);

    // Initialize a density matrix as a Hermitian matrix.
    GQCP::Orbital2DM<double> d {dim};

    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
            for (size_t k = 0; k < dim; k++) {
                for (size_t l = 0; l < dim; l++) {
                    d(i, j, k, l) = 1;
                }
            }
        }
    }

    // Initialize a reference value
    const double reference_expectation_value = 180.0;

    const double expectation_value = op.calculateExpectationValue(d);  // A scalar-StorageArray can be implicitly casted into its underlying scalar.
    BOOST_CHECK(std::abs(expectation_value - reference_expectation_value) < 1.0e-08);
}


/**
 *  Check if rotating with an identity transformation doesn't change the parameters.
 */
BOOST_AUTO_TEST_CASE(rotate_trivial) {

    const size_t dim = 2;

    // Set up toy two-electron integrals.
    auto op = toyTwoElectronIntegrals(dim);

    // Initialize an identity transformation.
    const auto U = GQCP::RTransformation<double>::Identity(dim);

    BOOST_CHECK(op.rotated(U).parameters().isApprox(op.parameters(), 1.0e-08));
}


/**
 *  Check the basis transformation formula using an other implementation (the old olsens code) from Ayers' Lab.
 */
BOOST_AUTO_TEST_CASE(transform_olsens) {

    // Set an example transformation and two-electron integrals.
    const size_t dim = 2;
    GQCP::SquareMatrix<double> T_matrix {dim};
    // clang-format off
    T_matrix << 1, 2,
                3, 4;
    // clang-format on
    const GQCP::RTransformation<double> T {T_matrix};


    auto g_par = GQCP::SquareRankFourTensor<double>::Zero(dim);
    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
            for (size_t k = 0; k < dim; k++) {
                for (size_t l = 0; l < dim; l++) {
                    g_par(i, j, k, l) = l + 2 * k + 4 * j + 8 * i;
                }
            }
        }
    }
    GQCP::ScalarRSQTwoElectronOperator<double> g {g_par};

    // Read in the reference and check the result.
    const auto g_transformed_par_ref = GQCP::SquareRankFourTensor<double>::FromFile("data/rotated_two_electron_integrals_olsens.data", dim);
    BOOST_CHECK(g.transformed(T).parameters().isApprox(g_transformed_par_ref, 1.0e-12));
}


/**
 *  Check whether or not the transformation method works as expected.
 */
BOOST_AUTO_TEST_CASE(transform) {

    const size_t dim = 2;

    // Initialize a test two-electron operator.
    auto op = toyTwoElectronIntegrals2(dim);

    // Initialize a test transformation.
    GQCP::SquareMatrix<double> T_matrix {dim};
    // clang-format off
    T_matrix << 2.0, 3.0,
                3.0, 4.0;
    // clang-format
    const GQCP::RTransformation<double> T {T_matrix};

    // Initialize the reference tensor.
    GQCP::SquareRankFourTensor<double> ref {dim};
    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
            for (size_t k = 0; k < dim; k++) {
                for (size_t l = 0; l < dim; l++) {
                    if ((i + j + k + l) == 0) {
                        ref(i, j, k, l) = 625.0;
                    }
                    if ((i + j + k + l) == 1) {
                        ref(i, j, k, l) = 875.0;
                    }
                    if ((i + j + k + l) == 2) {
                        ref(i, j, k, l) = 1225.0;
                    }
                    if ((i + j + k + l) == 3) {
                        ref(i, j, k, l) = 1715.0;
                    }
                    if ((i + j + k + l) == 4) {
                        ref(i, j, k, l) = 2401.0;
                    }
                }
            }
        }
    }

    op.transform(T);
    BOOST_CHECK(op.parameters().isApprox(ref, 1.0e-08));
}


/**
 * Check whether or not the jacobi rotation method works as expected
 */
BOOST_AUTO_TEST_CASE(transform_with_jacobi_matrix) {

    const size_t dim = 2;

    // Initialize a test two-electron operator.
    auto op = toyTwoElectronIntegrals2(dim);

    // Initialize the Jacobi rotation.
    GQCP::JacobiRotation J {1, 0, (boost::math::constants::pi<double>() / 2)};

    // Initialize the reference tensor.
    GQCP::SquareRankFourTensor<double> ref {dim};

    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
            for (size_t k = 0; k < dim; k++) {
                for (size_t l = 0; l < dim; l++) {
                    if ((i + j + k + l) % 2 == 0) {
                        ref(i, j, k, l) = 1.0;
                    }
                    if ((i + j + k + l) % 2 != 0) {
                        ref(i, j, k, l) = -1.0;
                    }
                }
            }
        }
    }

    op.rotate(J);
    BOOST_CHECK(op.parameters().isApprox(ref, 1.0e-08));
}


/**
 *  Check if antisymmetrizing two-electron integrals works as expected.
 * 
 *  The two-electron integrals under consideration are those from H2 in an STO-3G basisset.
 */
BOOST_AUTO_TEST_CASE(antisymmetrize) {

    // Prepare the two-electron repulsion integrals from the molecular Hamiltonian for H2.
    const auto molecule = GQCP::Molecule::HChain(2, 1.0);
    const GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> r_spinor_basis {molecule, "STO-3G"};
    const auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(r_spinor_basis, molecule);
    const auto& g = sq_hamiltonian.twoElectron();  // In chemist's notation.
    const auto g_A = g.antisymmetrized();

    // Antisymmetrize the chemist's two-electron integrals and check the results.
    const auto g_A_par = g_A.parameters();
    for (size_t p = 0; p < 2; p++) {
        for (size_t q = 0; q < 2; q++) {
            for (size_t r = 0; r < 2; r++) {
                for (size_t s = 0; s < 2; s++) {
                    BOOST_CHECK(std::abs(g_A_par(p, q, r, s) + g_A_par(r, q, p, s)) < 1.0e-12);
                    BOOST_CHECK(std::abs(g_A_par(p, q, r, s) + g_A_par(p, s, r, q)) < 1.0e-12);
                    BOOST_CHECK(std::abs(g_A_par(p, q, r, s) - g_A_par(r, s, p, q)) < 1.0e-12);
                }
            }
        }
    }


    // Convert the chemist's to physicist's two-electron integrals, antisymmetrize them and check the results.
    const auto V_A_par = g.convertedToPhysicistsNotation().antisymmetrized().parameters();
    for (size_t p = 0; p < 2; p++) {
        for (size_t q = 0; q < 2; q++) {
            for (size_t r = 0; r < 2; r++) {
                for (size_t s = 0; s < 2; s++) {
                    BOOST_CHECK(std::abs(V_A_par(p, q, r, s) + V_A_par(q, p, r, s)) < 1.0e-12);
                    BOOST_CHECK(std::abs(V_A_par(p, q, r, s) + V_A_par(p, q, s, r)) < 1.0e-12);
                    BOOST_CHECK(std::abs(V_A_par(p, q, r, s) - V_A_par(q, p, s, r)) < 1.0e-12);
                }
            }
        }
    }
}


/**
 *  Check if converting in-between chemist's and physicist's notation of two-electron integrals works as expected.
 * 
 *  The two-electron integrals under consideration are those from H2 in an STO-3G basisset.
 */
BOOST_AUTO_TEST_CASE(chemists_physicists) {

    // Prepare the two-electron repulsion integrals from the molecular Hamiltonian for H2.
    const auto molecule = GQCP::Molecule::HChain(2, 1.0);
    const GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> r_spinor_basis {molecule, "STO-3G"};
    const auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(r_spinor_basis, molecule);
    const auto& g = sq_hamiltonian.twoElectron();  // In chemist's notation.


    // Check if modifying chemist's integrals to chemist's integrals is a no-operation.
    const auto g_chemists = g.convertedToChemistsNotation();
    const auto& g_chemists_par = g_chemists.parameters();
    BOOST_CHECK(g_chemists.parameters().isApprox(g.parameters(), 1.0e-12));


    // Check if the conversion from chemist's to physicist's notation works as expected.
    const auto V_physicists = g_chemists.convertedToPhysicistsNotation();
    const auto& V_physicists_par = V_physicists.parameters();
    for (size_t p = 0; p < 2; p++) {
        for (size_t q = 0; q < 2; q++) {
            for (size_t r = 0; r < 2; r++) {
                for (size_t s = 0; s < 2; s++) {
                    BOOST_CHECK(std::abs(V_physicists_par(p, q, r, s) - g_chemists_par(p, r, q, s)) < 1.0e-12);
                }
            }
        }
    }


    // Check if modifying physicist's integrals to physicist's integrals is a no-operation.
    const auto V = V_physicists.convertedToPhysicistsNotation();
    BOOST_CHECK(V.parameters().isApprox(V_physicists.parameters(), 1.0e-12));


    // Check if the conversion from physicist's to chemist's notation works as expected.
    const auto g_again = V.convertedToChemistsNotation();
    BOOST_CHECK(g_again.parameters().isApprox(g.parameters(), 1.0e-12));
}
