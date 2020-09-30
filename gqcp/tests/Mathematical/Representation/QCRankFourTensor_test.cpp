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

#define BOOST_TEST_MODULE "QCRankFourTensor"

#include <boost/test/unit_test.hpp>

#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Basis/SpinorBasis/RSpinorBasis.hpp"
#include "Mathematical/Representation/QCRankFourTensor.hpp"
#include "Molecule/Molecule.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"

#include <boost/math/constants/constants.hpp>


/**
 *  Check the interface for QCRankFourTensor constructors.
 */
BOOST_AUTO_TEST_CASE(constructor) {

    GQCP::Tensor<double, 4> T1 {2, 2, 2, 2};
    T1.setZero();
    BOOST_CHECK_NO_THROW(GQCP::QCRankFourTensor<double> square_T1(T1));

    const GQCP::Tensor<double, 4> T2 {2, 1, 2, 2};
    BOOST_CHECK_THROW(GQCP::QCRankFourTensor<double> square_T2(T2), std::invalid_argument);  // not square
}


/**
 *  Check the basis transformation formula for a trivial case: T being a unit matrix.
 */
BOOST_AUTO_TEST_CASE(QCRankFourTensor_basisTransformInPlace_trivial) {

    const GQCP::TransformationMatrix<double> T = GQCP::TransformationMatrix<double>::Identity(3);

    GQCP::QCRankFourTensor<double> G {3};
    const auto G_copy = G;  // the reference
    G.basisTransform(T);

    BOOST_CHECK(G_copy.isApprox(G, 1.0e-12));
}


/**
 *  Check the basis transformation formula using an other implementation (the old olsens code) from Ayers' Lab.
 */
BOOST_AUTO_TEST_CASE(SQTwoElectronOperator_transform_olsens) {

    // Set an example transformation matrix and two-electron integrals tensor
    const size_t dim = 2;
    GQCP::TransformationMatrix<double> T {dim};
    // clang-format off
    T << 1, 2,
         3, 4;
    // clang-format on

    GQCP::QCRankFourTensor<double> G {dim};
    G.setZero();
    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
            for (size_t k = 0; k < dim; k++) {
                for (size_t l = 0; l < dim; l++) {
                    G(i, j, k, l) = l + 2 * k + 4 * j + 8 * i;
                }
            }
        }
    }
    G.basisTransform(T);


    // Read in the reference and check
    const GQCP::QCRankFourTensor<double> g_transformed_ref = GQCP::QCRankFourTensor<double>::FromFile("data/rotated_two_electron_integrals_olsens.data", dim);
    BOOST_CHECK(G.isApprox(g_transformed_ref, 1.0e-12));
}


/**
 *  Check if the code throws errors concerning non-unitary matrices being given to .rotate().
 */
BOOST_AUTO_TEST_CASE(QCRankFourTensor_rotate_throws) {

    // Create a random QCRankFourTensor
    const size_t dim = 3;
    GQCP::QCRankFourTensor<double> g {dim};
    g.setRandom();


    // Check if a non-unitary matrix as transformation matrix causes a throw
    const GQCP::TransformationMatrix<double> T = GQCP::TransformationMatrix<double>::Random(dim);  // chances are practically zero that a random matrix is unitary
    BOOST_CHECK_THROW(g.basisRotate(GQCP::TransformationMatrix<double>(T)), std::invalid_argument);


    // Check if a unitary matrix as transformation matrix is accepted
    GQCP::TransformationMatrix<double> U = GQCP::TransformationMatrix<double>::Identity(dim);
    g.basisRotate(U);
}


/**
 *  Check that rotating using JacobiRotationParameters is the same as using the corresponding Jacobi rotation matrix.
 */
BOOST_AUTO_TEST_CASE(QCRankFourTensor_basisRotateInPlace_JacobiRotationParameters) {

    // Create a random QCRankFourTensor
    const size_t dim = 5;
    GQCP::QCRankFourTensor<double> g {dim};
    g.setRandom();
    GQCP::QCRankFourTensor<double> G1 {g};
    GQCP::QCRankFourTensor<double> G2 {g};


    const GQCP::JacobiRotationParameters jacobi_rotation_parameters(4, 2, 56.81);
    const auto U = GQCP::TransformationMatrix<double>::FromJacobi(jacobi_rotation_parameters, dim);

    G1.basisRotate(jacobi_rotation_parameters);
    G2.basisRotate(U);


    BOOST_CHECK(G1.isApprox(G2, 1.0e-12));
}


/**
 *  Check if antisymmetrizing two-electron integrals works as expected.
 * 
 *  The two-electron integrals under consideration are those from H2 in an STO-3G basisset.
 */
BOOST_AUTO_TEST_CASE(antisymmetrize) {

    // Prepare the two-electron repulsion integrals from the molecular Hamiltonian for H2.
    const auto molecule = GQCP::Molecule::HChain(2, 1.0);
    const GQCP::RSpinorBasis<double, GQCP::GTOShell> r_spinor_basis {molecule, "STO-3G"};
    const auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(r_spinor_basis, molecule);
    const auto& g = sq_hamiltonian.twoElectron().parameters();  // in chemist's notation


    // Antisymmetrize the chemist's two-electron integrals and check the results.
    const auto g_A = g.antisymmetrized();
    for (size_t p = 0; p < 2; p++) {
        for (size_t q = 0; q < 2; q++) {
            for (size_t r = 0; r < 2; r++) {
                for (size_t s = 0; s < 2; s++) {
                    BOOST_CHECK(std::abs(g_A(p, q, r, s) + g_A(r, q, p, s)) < 1.0e-12);
                    BOOST_CHECK(std::abs(g_A(p, q, r, s) + g_A(p, s, r, q)) < 1.0e-12);
                    BOOST_CHECK(std::abs(g_A(p, q, r, s) - g_A(r, s, p, q)) < 1.0e-12);
                }
            }
        }
    }


    // Convert the chemist's to physicist's two-electron integrals, antisymmetrize them and check the results.
    const auto V_A = g.convertedToPhysicistsNotation().antisymmetrized();
    for (size_t p = 0; p < 2; p++) {
        for (size_t q = 0; q < 2; q++) {
            for (size_t r = 0; r < 2; r++) {
                for (size_t s = 0; s < 2; s++) {
                    BOOST_CHECK(std::abs(V_A(p, q, r, s) + V_A(q, p, r, s)) < 1.0e-12);
                    BOOST_CHECK(std::abs(V_A(p, q, r, s) + V_A(p, q, s, r)) < 1.0e-12);
                    BOOST_CHECK(std::abs(V_A(p, q, r, s) - V_A(q, p, s, r)) < 1.0e-12);
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
    const GQCP::RSpinorBasis<double, GQCP::GTOShell> r_spinor_basis {molecule, "STO-3G"};
    const auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(r_spinor_basis, molecule);
    const auto& g = sq_hamiltonian.twoElectron().parameters();  // in chemist's notation


    // Check if modifying chemist's integrals to chemist's integrals is a no-operation.
    const auto g_chemists = g.convertedToChemistsNotation();
    BOOST_CHECK(g_chemists.isApprox(g, 1.0e-12));


    // Check if the conversion from chemist's to physicist's notation works as expected.
    const auto V_physicists = g_chemists.convertedToPhysicistsNotation();
    for (size_t p = 0; p < 2; p++) {
        for (size_t q = 0; q < 2; q++) {
            for (size_t r = 0; r < 2; r++) {
                for (size_t s = 0; s < 2; s++) {
                    BOOST_CHECK(std::abs(V_physicists(p, q, r, s) - g_chemists(p, r, q, s)) < 1.0e-12);
                }
            }
        }
    }


    V_physicists.print();
    std::cout << std::endl
              << std::endl;

    g_chemists.print();


    // Check if modifying physicist's integrals to physicist's integrals is a no-operation.
    const auto V = V_physicists.convertedToPhysicistsNotation();
    BOOST_CHECK(V.isApprox(V_physicists, 1.0e-12));


    // Check if the conversion from physicist's to chemist's notation works as expected.
    const auto g_again = V.convertedToChemistsNotation();
    BOOST_CHECK(g_again.isApprox(g, 1.0e-12));
}
