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
#define BOOST_TEST_MODULE "AOBasis"

#include <boost/test/unit_test.hpp>

#include "Basis/AOBasis.hpp"

#include "Molecule/Molecule.hpp"


BOOST_AUTO_TEST_CASE ( AOBasis_constructor ) {

    // Check if we can construct an AOBasis object
    auto water = GQCP::Molecule::ReadXYZ("data/h2o.xyz");
    GQCP::AOBasis basis (water, "STO-3G");
}


BOOST_AUTO_TEST_CASE ( numberOfBasisFunctions ) {

    // Check the number of basis functions in water
    auto water = GQCP::Molecule::ReadXYZ("data/h2o.xyz");
    GQCP::AOBasis basis (water, "STO-3G");

    BOOST_CHECK_EQUAL(basis.numberOfBasisFunctions(), 7);
}


BOOST_AUTO_TEST_CASE( Szabo_integrals_h2_sto3g ) {

    // We will follow the example in Szabo, section 3.5.2, where it is stated that R = 1.4 a.u. = 0.740848 Angstrom
    auto h2 = GQCP::Molecule::ReadXYZ("data/h2_szabo.xyz");
    GQCP::AOBasis ao_basis (h2, "STO-3G");
    BOOST_CHECK_EQUAL(ao_basis.numberOfBasisFunctions(), 2);

    // Calculate some integrals
    auto S = ao_basis.calculateLibintOverlapIntegrals();
    auto T = ao_basis.calculateLibintKineticIntegrals();
    auto V = ao_basis.calculateLibintNuclearIntegrals();

    GQCP::OneElectronOperator<double> H_core = T + V;

    auto g = ao_basis.calculateLibintCoulombRepulsionIntegrals();

    // Fill in the reference values from Szabo
    GQCP::OneElectronOperator<double> ref_S (2);
    ref_S << 1.0,    0.6593,
             0.6593, 1.0;

    GQCP::OneElectronOperator<double> ref_T (2);
    ref_T << 0.7600, 0.2365,
             0.2365, 0.7600;

    GQCP::OneElectronOperator<double> ref_H_core (2);
    ref_H_core << -1.1204, -0.9584,
                  -0.9584, -1.1204;

    BOOST_CHECK(S.isApprox(ref_S, 1.0e-04));
    BOOST_CHECK(T.isApprox(ref_T, 1.0e-04));
    BOOST_CHECK(H_core.isApprox(ref_H_core, 1.0e-04));


    // The two-electron integrals in Szabo are given in chemist's notation, so this confirms that AO basis gives them in chemist's notation as well
    BOOST_CHECK(std::abs(g(0,0,0,0) - 0.7746) < 1.0e-04);
    BOOST_CHECK(std::abs(g(0,0,0,0) - g(1,1,1,1)) < 1.0e-12);

    BOOST_CHECK(std::abs(g(0,0,1,1) - 0.5697) < 1.0e-04);

    BOOST_CHECK(std::abs(g(1,0,0,0) - 0.4441) < 1.0e-04);
    BOOST_CHECK(std::abs(g(1,0,0,0) - g(1,1,1,0)) < 1.0e-12);

    BOOST_CHECK(std::abs(g(1,0,1,0) - 0.2970) < 1.0e-04);
}


BOOST_AUTO_TEST_CASE( HORTON_integrals_h2o_sto3g ) {

    // Set up an AO basis
    auto water = GQCP::Molecule::ReadXYZ("data/h2o.xyz");
    GQCP::AOBasis ao_basis (water, "STO-3G");
    auto nbf = ao_basis.numberOfBasisFunctions();


    // Calculate some integrals
    auto S = ao_basis.calculateLibintOverlapIntegrals();
    auto T = ao_basis.calculateLibintKineticIntegrals();
    auto V = ao_basis.calculateLibintNuclearIntegrals();

    auto g = ao_basis.calculateLibintCoulombRepulsionIntegrals();


    // Read in reference data from HORTON
    GQCP::OneElectronOperator<double> ref_S = GQCP::OneElectronOperator<double>::FromFile("data/h2o_sto-3g_overlap_horton.data", nbf, nbf);
    GQCP::OneElectronOperator<double> ref_T = GQCP::OneElectronOperator<double>::FromFile("data/h2o_sto-3g_kinetic_horton.data", nbf, nbf);
    GQCP::OneElectronOperator<double> ref_V = GQCP::OneElectronOperator<double>::FromFile("data/h2o_sto-3g_nuclear_horton.data", nbf, nbf);
    GQCP::TwoElectronOperator<double> ref_g = GQCP::TwoElectronOperator<double>::FromFile("data/h2o_sto-3g_coulomb_horton.data", nbf);


    // Check if the calculated integrals are close to those of HORTON
    BOOST_CHECK(S.isApprox(ref_S, 1.0e-07));
    BOOST_CHECK(T.isApprox(ref_T, 1.0e-07));
    BOOST_CHECK(V.isApprox(ref_V, 1.0e-07));
    BOOST_CHECK(g.isApprox(ref_g, 1.0e-06));
}


BOOST_AUTO_TEST_CASE ( libcint_vs_libint2_H2O_STO_3G ) {

    auto water = GQCP::Molecule::ReadXYZ("data/h2o.xyz");
    GQCP::AOBasis ao_basis (water, "STO-3G");

    const auto S_libcint = ao_basis.calculateLibcintOverlapIntegrals();
    const auto T_libcint = ao_basis.calculateLibcintKineticIntegrals();
    const auto V_libcint = ao_basis.calculateLibcintNuclearIntegrals();
    const auto dipole_libcint = ao_basis.calculateLibcintDipoleIntegrals();
    const auto g_libcint = ao_basis.calculateLibcintCoulombRepulsionIntegrals();

    const auto S_libint2 = ao_basis.calculateLibintOverlapIntegrals();
    const auto T_libint2 = ao_basis.calculateLibintKineticIntegrals();
    const auto V_libint2 = ao_basis.calculateLibintNuclearIntegrals();
    const auto dipole_libint2 = ao_basis.calculateLibintDipoleIntegrals();
    const auto g_libint2 = ao_basis.calculateLibintCoulombRepulsionIntegrals();
    
    BOOST_CHECK(S_libcint.isApprox(S_libint2, 1.0e-08));
    BOOST_CHECK(T_libcint.isApprox(T_libint2, 1.0e-08));
    BOOST_CHECK(V_libcint.isApprox(V_libint2, 1.0e-08));
    for (size_t i = 0; i < 3; i++) {
        BOOST_CHECK(dipole_libcint[i].isApprox(dipole_libint2[i], 1.0e-08));
    }
    BOOST_CHECK(g_libcint.isApprox(g_libint2, 1.0e-08));

}


BOOST_AUTO_TEST_CASE ( libcint_vs_libint2_dipole_origin ) {

    auto water = GQCP::Molecule::ReadXYZ("data/h2o.xyz");
    GQCP::AOBasis ao_basis (water, "STO-3G");

    GQCP::Vector<double, 3> origin;
    origin << 0.0, 1.0, -0.5;

    const auto dipole_libcint = ao_basis.calculateLibcintDipoleIntegrals(origin);
    const auto dipole_libint2 = ao_basis.calculateLibintDipoleIntegrals(origin);

    for (size_t i = 0; i < 3; i++) {
        BOOST_CHECK(dipole_libcint[i].isApprox(dipole_libint2[i], 1.0e-08));
    }
}
