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
#define BOOST_TEST_MODULE "IntegralCalculator"

#include <boost/test/unit_test.hpp>

#include "Basis/Integrals/IntegralCalculator.hpp"
#include "Basis/ScalarBasis/ScalarBasis.hpp"
#include "Molecule/Molecule.hpp"
#include "Operator/FirstQuantized/Operator.hpp"


/**
 *  Check integrals calculated by Libint with reference values in Szabo.
 */
BOOST_AUTO_TEST_CASE( Szabo_integrals_h2_sto3g ) {

    // In Szabo, section 3.5.2, we read that the internuclear distance R = 1.4 a.u. = 0.740848 Angstrom
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2_szabo.xyz");
    const GQCP::ScalarBasis<GQCP::GTOShell> scalar_basis (molecule, "STO-3G");
    BOOST_CHECK_EQUAL(scalar_basis.numberOfBasisFunctions(), 2);


    // Let Libint2 calculate some integrals
    const auto S = GQCP::IntegralCalculator::calculateLibintIntegrals(GQCP::Operator::Overlap(), scalar_basis);
    const auto T = GQCP::IntegralCalculator::calculateLibintIntegrals(GQCP::Operator::Kinetic(), scalar_basis);
    const auto V = GQCP::IntegralCalculator::calculateLibintIntegrals(GQCP::Operator::NuclearAttraction(molecule), scalar_basis);
    const GQCP::QCMatrix<double> H_core = T + V;

    const auto g = GQCP::IntegralCalculator::calculateLibintIntegrals(GQCP::Operator::Coulomb(), scalar_basis);


    // Check the one-electron integrals with the reference
    GQCP::QCMatrix<double> ref_S (2);
    ref_S << 1.0,    0.6593,
             0.6593, 1.0;

    GQCP::QCMatrix<double> ref_T (2);
    ref_T << 0.7600, 0.2365,
             0.2365, 0.7600;

    GQCP::QCMatrix<double> ref_H_core (2);
    ref_H_core << -1.1204, -0.9584,
                  -0.9584, -1.1204;

    BOOST_CHECK(S.isApprox(ref_S, 1.0e-04));
    BOOST_CHECK(T.isApprox(ref_T, 1.0e-04));
    BOOST_CHECK(H_core.isApprox(ref_H_core, 1.0e-04));


    // Check the two-electron integrals with the reference
    // The two-electron integrals in Szabo are given in chemist's notation, so this confirms that AO basis gives them in chemist's notation as well
    BOOST_CHECK(std::abs(g(0,0,0,0) - 0.7746) < 1.0e-04);
    BOOST_CHECK(std::abs(g(0,0,0,0) - g(1,1,1,1)) < 1.0e-12);

    BOOST_CHECK(std::abs(g(0,0,1,1) - 0.5697) < 1.0e-04);

    BOOST_CHECK(std::abs(g(1,0,0,0) - 0.4441) < 1.0e-04);
    BOOST_CHECK(std::abs(g(1,0,0,0) - g(1,1,1,0)) < 1.0e-12);

    BOOST_CHECK(std::abs(g(1,0,1,0) - 0.2970) < 1.0e-04);
}


/**
 *  Check integrals calculated by Libint with reference values from HORTON.
 */
BOOST_AUTO_TEST_CASE( HORTON_integrals_h2o_sto3g ) {

    // Set up an AO basis
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o.xyz");
    const GQCP::ScalarBasis<GQCP::GTOShell> scalar_basis (molecule, "STO-3G");
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
BOOST_AUTO_TEST_CASE ( libcint_vs_libint2_H2O_STO_3G ) {

    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o.xyz");
    const GQCP::ScalarBasis<GQCP::GTOShell> scalar_basis (molecule, "STO-3G");


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
BOOST_AUTO_TEST_CASE ( libcint_vs_libint2_dipole_origin ) {

    auto molecule = GQCP::Molecule::ReadXYZ("data/h2o.xyz");
    GQCP::ScalarBasis<GQCP::GTOShell> scalar_basis (molecule, "STO-3G");

    GQCP::Vector<double, 3> origin;
    origin << 0.0, 1.0, -0.5;

    const auto dipole_libint2 = GQCP::IntegralCalculator::calculateLibintIntegrals(GQCP::Operator::ElectronicDipole(origin), scalar_basis);
    const auto dipole_libcint = GQCP::IntegralCalculator::calculateLibcintIntegrals(GQCP::Operator::ElectronicDipole(origin), scalar_basis);

    for (size_t i = 0; i < 3; i++) {
        BOOST_CHECK(dipole_libcint[i].isApprox(dipole_libint2[i], 1.0e-08));
    }
}



