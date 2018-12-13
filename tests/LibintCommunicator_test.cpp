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
#define BOOST_TEST_MODULE "LibintCommunicator"


#include "LibintCommunicator.hpp"

#include "utilities/io.hpp"
#include "utilities/linalg.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain



BOOST_AUTO_TEST_CASE ( atoms_interface ) {

    std::vector<GQCP::Atom> GQCP_atoms = {
        {1, 0, 3, 0},
        {2, 0, 0, 4},
        {3, 3, 0, 0},
        {4, 0, 0, 5}
    };

    std::vector<libint2::Atom> ref_libint_atoms = {
        {1, 0, 3, 0},
        {2, 0, 0, 4},
        {3, 3, 0, 0},
        {4, 0, 0, 5}
    };


    // Use the Libint interface to obtain a std::vector<libint2::Atom> from the GQCP ones
    auto test_libint_atoms = GQCP::LibintCommunicator::get().interface(GQCP_atoms);


    /**
     *  Implement a function object to compare (libint_atom) == (libint_atom)
     */
    struct LibintAtomEqual {
        double tolerance;
        explicit LibintAtomEqual(double tolerance) : tolerance (tolerance) {};

        bool operator()(const libint2::Atom& lhs, const libint2::Atom& rhs) {
            return (lhs.atomic_number == rhs.atomic_number) &&
                   (std::abs(lhs.x - rhs.x) < tolerance) &&
                   (std::abs(lhs.y - rhs.y) < tolerance) &&
                   (std::abs(lhs.z - rhs.z) < tolerance);
        }
    };


    // Check if the interfacing between the Atom types works
    BOOST_CHECK((ref_libint_atoms.size() == test_libint_atoms.size()) &&
                std::equal(ref_libint_atoms.begin(), ref_libint_atoms.end(), test_libint_atoms.begin(), LibintAtomEqual(1.0e-08)));
}


BOOST_AUTO_TEST_CASE( Szabo_integrals_h2_sto3g ) {

    // We will follow the example in Szabo, section 3.5.2, where it is stated that R = 1.4 a.u. = 0.740848 Angstrom
    auto h2 = GQCP::Molecule::Readxyz("../tests/data/h2_szabo.xyz");
    GQCP::AOBasis basis (h2, "STO-3G");
    BOOST_CHECK_EQUAL(basis.get_number_of_basis_functions(), 2);

    // Calculate some integrals
    auto S = GQCP::LibintCommunicator::get().calculateOverlapIntegrals(basis);
    auto T = GQCP::LibintCommunicator::get().calculateKineticIntegrals(basis);
    auto V = GQCP::LibintCommunicator::get().calculateNuclearIntegrals(basis);

    auto H_core = GQCP::OneElectronOperator(T.get_matrix_representation() + V.get_matrix_representation());

    auto g = GQCP::LibintCommunicator::get().calculateCoulombRepulsionIntegrals(basis);


    // Fill in the reference values from Szabo
    Eigen::MatrixXd ref_S (2, 2);
    ref_S << 1.0,    0.6593,
             0.6593, 1.0;

    Eigen::MatrixXd ref_T (2, 2);
    ref_T << 0.7600, 0.2365,
             0.2365, 0.7600;

    Eigen::MatrixXd ref_H_core (2, 2);
    ref_H_core << -1.1204, -0.9584,
                  -0.9584, -1.1204;

    BOOST_CHECK(S.get_matrix_representation().isApprox(ref_S, 1.0e-04));
    BOOST_CHECK(T.get_matrix_representation().isApprox(ref_T, 1.0e-04));
    BOOST_CHECK(H_core.get_matrix_representation().isApprox(ref_H_core, 1.0e-04));


    // The two-electron integrals in Szabo are given in chemist's notation, so this confirms that the LibintCommunicator gives them in chemist's notation as well
    BOOST_CHECK(std::abs(g.get_matrix_representation()(0,0,0,0) - 0.7746) < 1.0e-04);
    BOOST_CHECK(std::abs(g.get_matrix_representation()(0,0,0,0) - g.get_matrix_representation()(1,1,1,1)) < 1.0e-12);

    BOOST_CHECK(std::abs(g.get_matrix_representation()(0,0,1,1) - 0.5697) < 1.0e-04);

    BOOST_CHECK(std::abs(g.get_matrix_representation()(1,0,0,0) - 0.4441) < 1.0e-04);
    BOOST_CHECK(std::abs(g.get_matrix_representation()(1,0,0,0) - g.get_matrix_representation()(1,1,1,0)) < 1.0e-12);

    BOOST_CHECK(std::abs(g.get_matrix_representation()(1,0,1,0) - 0.2970) < 1.0e-04);
}


BOOST_AUTO_TEST_CASE( HORTON_integrals_h2o_sto3g ) {

    // Set up a basis
    auto water = GQCP::Molecule::Readxyz("../tests/data/h2o.xyz");
    GQCP::AOBasis basis (water, "STO-3G");
    auto nbf = basis.get_number_of_basis_functions();


    // Calculate some integrals
    libint2::BasisSet basisset = basis.get_basis_functions();
    auto S = GQCP::LibintCommunicator::get().calculateOverlapIntegrals(basis);
    auto T = GQCP::LibintCommunicator::get().calculateKineticIntegrals(basis);
    auto V = GQCP::LibintCommunicator::get().calculateNuclearIntegrals(basis);

    auto g = GQCP::LibintCommunicator::get().calculateCoulombRepulsionIntegrals(basis);


    // Read in reference data from HORTON
    Eigen::MatrixXd ref_S (nbf, nbf);
    Eigen::MatrixXd ref_T (nbf, nbf);
    Eigen::MatrixXd ref_V (nbf, nbf);
    Eigen::Tensor<double, 4> ref_g (nbf, nbf, nbf, nbf);

    GQCP::readArrayFromFile("../tests/data/h2o_sto-3g_overlap_horton.data", ref_S);
    GQCP::readArrayFromFile("../tests/data/h2o_sto-3g_kinetic_horton.data", ref_T);
    GQCP::readArrayFromFile("../tests/data/h2o_sto-3g_nuclear_horton.data", ref_V);
    GQCP::readArrayFromFile("../tests/data/h2o_sto-3g_coulomb_horton.data", ref_g);


    // Check if the calculated integrals are equal to those of HORTON
    BOOST_CHECK(S.get_matrix_representation().isApprox(ref_S, 1.0e-08));
    BOOST_CHECK(T.get_matrix_representation().isApprox(ref_T, 1.0e-08));
    BOOST_CHECK(V.get_matrix_representation().isApprox(ref_V, 1.0e-08));
    BOOST_CHECK(GQCP::areEqual(g.get_matrix_representation(), ref_g, 1.0e-06));
}
