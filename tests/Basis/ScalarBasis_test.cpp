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
#define BOOST_TEST_MODULE "ScalarBasis"

#include <boost/test/unit_test.hpp>

#include "Basis/ScalarBasis.hpp"
#include "Molecule/Molecule.hpp"


BOOST_AUTO_TEST_CASE ( AOBasis_constructor ) {

    // Check if we can construct an ScalarBasis<GTOShell> object
    auto water = GQCP::Molecule::ReadXYZ("data/h2o.xyz");
    GQCP::ScalarBasis<GQCP::GTOShell> basis (water, "STO-3G");
}


BOOST_AUTO_TEST_CASE ( numberOfBasisFunctions ) {

    // Check the number of basis functions in water
    auto water = GQCP::Molecule::ReadXYZ("data/h2o.xyz");
    GQCP::ScalarBasis<GQCP::GTOShell> basis (water, "STO-3G");

    BOOST_CHECK_EQUAL(basis.numberOfBasisFunctions(), 7);
}


/**
 *  Check integrals calculated by Libint with reference values in Szabo
 */
BOOST_AUTO_TEST_CASE( Szabo_integrals_h2_sto3g ) {

    // In Szabo, section 3.5.2, we read that the internuclear distance R = 1.4 a.u. = 0.740848 Angstrom
    const auto h2 = GQCP::Molecule::ReadXYZ("data/h2_szabo.xyz");
    const GQCP::ScalarBasis<GQCP::GTOShell> ao_basis (h2, "STO-3G");
    BOOST_CHECK_EQUAL(ao_basis.numberOfBasisFunctions(), 2);


    // Let Libint2 calculate some integrals
    const auto S = ao_basis.calculateLibintOverlapIntegrals();
    const auto T = ao_basis.calculateLibintKineticIntegrals();
    const auto V = ao_basis.calculateLibintNuclearIntegrals();

    const GQCP::ChemicalMatrix<double> H_core = T + V;

    const auto g = ao_basis.calculateLibintCoulombRepulsionIntegrals();


    // Check the one-electron integrals with the reference
    GQCP::ChemicalMatrix<double> ref_S (2);
    ref_S << 1.0,    0.6593,
             0.6593, 1.0;

    GQCP::ChemicalMatrix<double> ref_T (2);
    ref_T << 0.7600, 0.2365,
             0.2365, 0.7600;

    GQCP::ChemicalMatrix<double> ref_H_core (2);
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


BOOST_AUTO_TEST_CASE( HORTON_integrals_h2o_sto3g ) {

    // Set up an AO basis
    const auto water = GQCP::Molecule::ReadXYZ("data/h2o.xyz");
    const GQCP::ScalarBasis<GQCP::GTOShell> ao_basis (water, "STO-3G");
    const auto nbf = ao_basis.numberOfBasisFunctions();


    // Calculate some integrals
    const auto S = ao_basis.calculateLibintOverlapIntegrals();
    const auto T = ao_basis.calculateLibintKineticIntegrals();
    const auto V = ao_basis.calculateLibintNuclearIntegrals();

    const auto g = ao_basis.calculateLibintCoulombRepulsionIntegrals();


    // Read in reference data from HORTON
    const auto ref_S = GQCP::ChemicalMatrix<double>::FromFile("data/h2o_sto-3g_overlap_horton.data", nbf, nbf);
    const auto ref_T = GQCP::ChemicalMatrix<double>::FromFile("data/h2o_sto-3g_kinetic_horton.data", nbf, nbf);
    const auto ref_V = GQCP::ChemicalMatrix<double>::FromFile("data/h2o_sto-3g_nuclear_horton.data", nbf, nbf);
    const auto ref_g = GQCP::ChemicalRankFourTensor<double>::FromFile("data/h2o_sto-3g_coulomb_horton.data", nbf);


    // Check if the calculated integrals are close to those of HORTON
    BOOST_CHECK(S.isApprox(ref_S, 1.0e-07));
    BOOST_CHECK(T.isApprox(ref_T, 1.0e-07));
    BOOST_CHECK(V.isApprox(ref_V, 1.0e-07));
    BOOST_CHECK(g.isApprox(ref_g, 1.0e-06));
}


BOOST_AUTO_TEST_CASE ( libcint_vs_libint2_H2O_STO_3G ) {

    auto water = GQCP::Molecule::ReadXYZ("data/h2o.xyz");
    GQCP::ScalarBasis<GQCP::GTOShell> ao_basis (water, "STO-3G");

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
    GQCP::ScalarBasis<GQCP::GTOShell> ao_basis (water, "STO-3G");

    GQCP::Vector<double, 3> origin;
    origin << 0.0, 1.0, -0.5;

    const auto dipole_libcint = ao_basis.calculateLibcintDipoleIntegrals(origin);
    const auto dipole_libint2 = ao_basis.calculateLibintDipoleIntegrals(origin);

    for (size_t i = 0; i < 3; i++) {
        BOOST_CHECK(dipole_libcint[i].isApprox(dipole_libint2[i], 1.0e-08));
    }
}


BOOST_AUTO_TEST_CASE ( dissociatedMoleculeBasis ) {

    // Test if we can succesfully initialize NO+ at long intra molecular distance
    auto N = GQCP::Nucleus(7, 3.5, 0, 0);
    auto O = GQCP::Nucleus(8, -3.5, 0, 0);
    std::vector<GQCP::Nucleus> nuclei {N,O};
    auto NO = GQCP::Molecule(nuclei, +1);
    GQCP::ScalarBasis<GQCP::GTOShell> basis (NO, "STO-3G");

    BOOST_CHECK_NO_THROW(GQCP::ScalarBasis<GQCP::GTOShell> basis (NO, "STO-3G"));
}
