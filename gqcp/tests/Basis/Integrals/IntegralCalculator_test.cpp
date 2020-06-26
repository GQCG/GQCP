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

#define BOOST_TEST_MODULE "IntegralCalculator"

#include <boost/test/unit_test.hpp>

#include "Basis/Integrals/IntegralCalculator.hpp"
#include "Basis/ScalarBasis/ScalarBasis.hpp"
#include "Molecule/Molecule.hpp"
#include "Operator/FirstQuantized/Operator.hpp"


/**
 *  Check integrals calculated by Libint with reference values in Szabo.
 */
BOOST_AUTO_TEST_CASE(Szabo_integrals_h2_sto3g) {

    // In Szabo, section 3.5.2, we read that the internuclear distance R = 1.4 a.u. = 0.740848 Angstrom
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2_szabo.xyz");
    const GQCP::ScalarBasis<GQCP::GTOShell> scalar_basis {molecule, "STO-3G"};
    BOOST_CHECK_EQUAL(scalar_basis.numberOfBasisFunctions(), 2);


    // Let Libint2 calculate some integrals
    const auto S = GQCP::IntegralCalculator::calculateLibintIntegrals(GQCP::Operator::Overlap(), scalar_basis);
    const auto T = GQCP::IntegralCalculator::calculateLibintIntegrals(GQCP::Operator::Kinetic(), scalar_basis);
    const auto V = GQCP::IntegralCalculator::calculateLibintIntegrals(GQCP::Operator::NuclearAttraction(molecule), scalar_basis);
    const GQCP::QCMatrix<double> H_core = T + V;

    const auto g = GQCP::IntegralCalculator::calculateLibintIntegrals(GQCP::Operator::Coulomb(), scalar_basis);


    // Check the one-electron integrals with the reference
    GQCP::QCMatrix<double> ref_S {2};
    // clang-format off
    ref_S << 1.0,    0.6593,
             0.6593, 1.0;
    // clang-format on

    GQCP::QCMatrix<double> ref_T {2};
    // clang-format off
    ref_T << 0.7600, 0.2365,
             0.2365, 0.7600;
    // clang-format on

    GQCP::QCMatrix<double> ref_H_core {2};
    // clang-format off
    ref_H_core << -1.1204, -0.9584,
                  -0.9584, -1.1204;
    // clang-format on

    BOOST_CHECK(S.isApprox(ref_S, 1.0e-04));
    BOOST_CHECK(T.isApprox(ref_T, 1.0e-04));
    BOOST_CHECK(H_core.isApprox(ref_H_core, 1.0e-04));


    // Check the two-electron integrals with the reference
    // The two-electron integrals in Szabo are given in chemist's notation, so this confirms that AO basis gives them in chemist's notation as well
    BOOST_CHECK(std::abs(g(0, 0, 0, 0) - 0.7746) < 1.0e-04);
    BOOST_CHECK(std::abs(g(0, 0, 0, 0) - g(1, 1, 1, 1)) < 1.0e-12);

    BOOST_CHECK(std::abs(g(0, 0, 1, 1) - 0.5697) < 1.0e-04);

    BOOST_CHECK(std::abs(g(1, 0, 0, 0) - 0.4441) < 1.0e-04);
    BOOST_CHECK(std::abs(g(1, 0, 0, 0) - g(1, 1, 1, 0)) < 1.0e-12);

    BOOST_CHECK(std::abs(g(1, 0, 1, 0) - 0.2970) < 1.0e-04);
}


/**
 *  Check integrals calculated by Libint with reference values from HORTON.
 */
BOOST_AUTO_TEST_CASE(HORTON_integrals_h2o_sto3g) {

    // Set up an AO basis
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o.xyz");
    const GQCP::ScalarBasis<GQCP::GTOShell> scalar_basis {molecule, "STO-3G"};
    const auto nbf = scalar_basis.numberOfBasisFunctions();


    // Calculate some integrals
    const auto S = GQCP::IntegralCalculator::calculateLibintIntegrals(GQCP::Operator::Overlap(), scalar_basis);
    const auto T = GQCP::IntegralCalculator::calculateLibintIntegrals(GQCP::Operator::Kinetic(), scalar_basis);
    const auto V = GQCP::IntegralCalculator::calculateLibintIntegrals(GQCP::Operator::NuclearAttraction(molecule), scalar_basis);
    const auto g = GQCP::IntegralCalculator::calculateLibintIntegrals(GQCP::Operator::Coulomb(), scalar_basis);


    // Read in reference data from HORTON
    const auto ref_S = GQCP::QCMatrix<double>::FromFile("data/h2o_sto-3g_overlap_horton.data", nbf, nbf);
    const auto ref_T = GQCP::QCMatrix<double>::FromFile("data/h2o_sto-3g_kinetic_horton.data", nbf, nbf);
    const auto ref_V = GQCP::QCMatrix<double>::FromFile("data/h2o_sto-3g_nuclear_horton.data", nbf, nbf);
    const auto ref_g = GQCP::QCRankFourTensor<double>::FromFile("data/h2o_sto-3g_coulomb_horton.data", nbf);


    // Check if the calculated integrals are close to those of HORTON
    BOOST_CHECK(S.isApprox(ref_S, 1.0e-07));
    BOOST_CHECK(T.isApprox(ref_T, 1.0e-07));
    BOOST_CHECK(V.isApprox(ref_V, 1.0e-07));
    BOOST_CHECK(g.isApprox(ref_g, 1.0e-06));
}


/**
 *  Check the calculation of some integrals between Libint2 and libcint.
 */
BOOST_AUTO_TEST_CASE(libcint_vs_libint2_H2O_STO_3G) {

    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o.xyz");
    const GQCP::ScalarBasis<GQCP::GTOShell> scalar_basis {molecule, "STO-3G"};


    const auto S_libint2 = GQCP::IntegralCalculator::calculateLibintIntegrals(GQCP::Operator::Overlap(), scalar_basis);
    const auto T_libint2 = GQCP::IntegralCalculator::calculateLibintIntegrals(GQCP::Operator::Kinetic(), scalar_basis);
    const auto V_libint2 = GQCP::IntegralCalculator::calculateLibintIntegrals(GQCP::Operator::NuclearAttraction(molecule), scalar_basis);
    const auto dipole_libint2 = GQCP::IntegralCalculator::calculateLibintIntegrals(GQCP::Operator::ElectronicDipole(), scalar_basis);
    const auto g_libint2 = GQCP::IntegralCalculator::calculateLibintIntegrals(GQCP::Operator::Coulomb(), scalar_basis);

    const auto S_libcint = GQCP::IntegralCalculator::calculateLibcintIntegrals(GQCP::Operator::Overlap(), scalar_basis);
    const auto T_libcint = GQCP::IntegralCalculator::calculateLibcintIntegrals(GQCP::Operator::Kinetic(), scalar_basis);
    const auto V_libcint = GQCP::IntegralCalculator::calculateLibcintIntegrals(GQCP::Operator::NuclearAttraction(molecule), scalar_basis);
    const auto dipole_libcint = GQCP::IntegralCalculator::calculateLibcintIntegrals(GQCP::Operator::ElectronicDipole(), scalar_basis);
    const auto g_libcint = GQCP::IntegralCalculator::calculateLibcintIntegrals(GQCP::Operator::Coulomb(), scalar_basis);


    BOOST_CHECK(S_libcint.isApprox(S_libint2, 1.0e-08));
    BOOST_CHECK(T_libcint.isApprox(T_libint2, 1.0e-08));
    BOOST_CHECK(V_libcint.isApprox(V_libint2, 1.0e-08));
    for (size_t i = 0; i < 3; i++) {
        BOOST_CHECK(dipole_libcint[i].isApprox(dipole_libint2[i], 1.0e-08));
    }
    BOOST_CHECK(g_libcint.isApprox(g_libint2, 1.0e-08));
}


/**
 *  Check the dipole integrals between libcint and libint2 for an origin different from zero.
 */
BOOST_AUTO_TEST_CASE(libcint_vs_libint2_dipole_origin) {

    auto molecule = GQCP::Molecule::ReadXYZ("data/h2o.xyz");
    GQCP::ScalarBasis<GQCP::GTOShell> scalar_basis {molecule, "STO-3G"};

    GQCP::Vector<double, 3> origin;
    origin << 0.0, 1.0, -0.5;

    const auto dipole_libint2 = GQCP::IntegralCalculator::calculateLibintIntegrals(GQCP::Operator::ElectronicDipole(origin), scalar_basis);
    const auto dipole_libcint = GQCP::IntegralCalculator::calculateLibcintIntegrals(GQCP::Operator::ElectronicDipole(origin), scalar_basis);

    for (size_t i = 0; i < 3; i++) {
        BOOST_CHECK(dipole_libcint[i].isApprox(dipole_libint2[i], 1.0e-08));
    }
}


/**
 *  Check for the correctness of the angular momentum integrals that libcint produces.
 */
// BOOST_AUTO_TEST_CASE(libcint_angular_momentum) {

//     // Set up the example scalar basis.
//     const auto molecule = GQCP::Molecule::ReadXYZ("data/h2_szabo.xyz");
//     const GQCP::ScalarBasis<GQCP::GTOShell> scalar_basis {molecule, "STO-3G"};

//     // const auto s = scalar_basis.basisFunctions()[0];

//     // std::cout << "1s_L: " << std::endl
//     //           << s.description() << std::endl
//     //           << std::endl;


//     const auto L = GQCP::IntegralCalculator::calculateLibcintIntegrals(GQCP::Operator::AngularMomentum(), scalar_basis);

//     std::cout << "Lx integrals: " << std::endl
//               << L[0] << std::endl
//               << std::endl;
// }


#include <boost/math/constants/constants.hpp>
#include "Mathematical/CartesianDirection.hpp"


/**
 *  Calculate the McMurchie-Davidson expansion coefficient E_t^{ij} for one of the Cartesian directions.
 * 
 *  K: K_x
 *  L: L_x
 * 
 *  a,b exponents
 */
double calculateE(const size_t i, const size_t j, const int t, const double K, const double L, const double a, const double b) {

    // std::cout << "Calculating E^{ij}_t" << i << ' ' << j << ' ' << t << std::endl;

    const double Delta = K - L;

    // Check if t is out of bounds: 0 <= t <= i+j.
    if ((t < 0) || (t > (i + j))) {
        return 0.0;
    }


    // Provide the base recurrence case.
    else if ((t == 0) && (i == 0) && (j == 0)) {
        const double mu = a * b / (a + b);

        return std::exp(-mu * std::pow(Delta, 2));
    }

    // Do the recurrence for E^{i+1, j}_t.
    else if ((j == 0)) {
        // Do the actual recursion, but calculate the equal terms in the two cases first.
        const double p = a + b;

        double value = 1.0 / (2 * p) * calculateE(i - 1, j, t - 1, K, L, a, b) + (t + 1) * calculateE(i - 1, j, t + 1, K, L, a, b);

        return value -= b / p * Delta * calculateE(i - 1, j, t, K, L, a, b);
    }

    else {  // do the recurrence for E^{i, j+1}_t

        // Do the actual recursion, but calculate the equal terms in the two cases first.
        const double p = a + b;

        double value = 1.0 / (2 * p) * calculateE(i, j - 1, t - 1, K, L, a, b) + (t + 1) * calculateE(i, j - 1, t + 1, K, L, a, b);

        return value += a / p * Delta * calculateE(i, j - 1, t, K, L, a, b);
    }
}


BOOST_AUTO_TEST_CASE(overlap) {

    // Set up an AO basis
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o.xyz");
    const GQCP::ScalarBasis<GQCP::GTOShell> scalar_basis {molecule, "STO-3G"};
    const auto S_libint2 = GQCP::IntegralCalculator::calculateLibintIntegrals(GQCP::Operator::Overlap(), scalar_basis);

    std::cout << "S: " << std::endl
              << S_libint2 << std::endl
              << std::endl;

    std::vector<double> integrals;
    for (const auto sh1 : scalar_basis.shellSet().asVector()) {
        const auto K = sh1.nucleus().position();

        for (const auto sh2 : scalar_basis.shellSet().asVector()) {
            const auto L = sh2.nucleus().position();

            for (const auto cartesian_exponents1 : sh1.generateCartesianExponents()) {

                for (const auto cartesian_exponents2 : sh2.generateCartesianExponents()) {
                    std::cout << cartesian_exponents1.description() << ' ' << cartesian_exponents2.description() << std::endl;

                    double integral = 0.0;

                    for (size_t c1 = 0; c1 < sh1.contractionSize(); c1++) {
                        for (size_t c2 = 0; c2 < sh2.contractionSize(); c2++) {
                            const auto a = sh1.gaussianExponents()[c1];
                            const auto b = sh2.gaussianExponents()[c2];

                            const auto d1 = sh1.contractionCoefficients()[c1];
                            const auto d2 = sh2.contractionCoefficients()[c2];

                            const auto p = a + b;
                            double primitive_integral = std::pow(boost::math::constants::pi<double>() / p, 1.5);

                            for (const auto& direction : {GQCP::CartesianDirection::x, GQCP::CartesianDirection::y, GQCP::CartesianDirection::z}) {
                                primitive_integral *= calculateE(cartesian_exponents1.value(direction), cartesian_exponents2.value(direction), 0, K(direction), L(direction), a, b);
                            }

                            integral += d1 * d2 * primitive_integral;
                        }
                    }

                    integrals.push_back(integral);
                }
            }
        }
    }


    for (const auto integral : integrals) {
        std::cout << integral << std::endl;
    }

    std::cout << integrals.size() << std::endl;
}
