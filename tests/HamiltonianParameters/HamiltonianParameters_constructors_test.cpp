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
#define BOOST_TEST_MODULE "HamiltonianParameters_constructors"


#include "HamiltonianParameters/HamiltonianParameters_constructors.hpp"

#include <cpputil.hpp>

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain



BOOST_AUTO_TEST_CASE ( constructMolecularHamiltonianParameters ) {

    // Set up a basis
    GQCP::Molecule water ("../tests/data/h2o.xyz");
    auto ao_basis_sptr = std::make_shared<GQCP::AOBasis>(water, "STO-3G");

    
    // Check if we can construct the molecular Hamiltonian parameters
    auto mol_ham_par = GQCP::constructMolecularHamiltonianParameters(ao_basis_sptr);
}


BOOST_AUTO_TEST_CASE ( constructRandomHamiltonianParameters ) {

    auto random_ham_par = GQCP::constructRandomHamiltonianParameters(5);
}


BOOST_AUTO_TEST_CASE ( FCIDUMP_reader ) {
    
    auto fcidump_ham_par = GQCP::readFCIDUMPFile("../tests/data/beh_cation_631g_caitlin.FCIDUMP");
    
    // Check if the one-electron integrals are read in correctly from a previous implementation
    GQCP::OneElectronOperator h_SO = fcidump_ham_par.get_h();
    
    BOOST_CHECK(std::abs(h_SO(0,0) - (-8.34082)) < 1.0e-5);
    BOOST_CHECK(std::abs(h_SO(5,1) - 0.381418) < 1.0e-6);
    BOOST_CHECK(std::abs(h_SO(14,0) - 0.163205) < 1.0e-6);
    BOOST_CHECK(std::abs(h_SO(13,6) - (-5.53204e-16)) < 1.0e-16);
    BOOST_CHECK(std::abs(h_SO(15,11) - (-0.110721)) < 1.0e-6);
    
    
    // Check if the two-electron integrals are read in correctly from a previous implementation
    GQCP::TwoElectronOperator g_SO = fcidump_ham_par.get_g();
    
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
    auto fcidump_ham_par = GQCP::readFCIDUMPFile("../tests/data/h2_psi4_horton.FCIDUMP");
    
    GQCP::TwoElectronOperator g_SO = fcidump_ham_par.get_g();
    BOOST_CHECK(std::abs(g_SO(6,5,1,0) - 0.0533584656) <  1.0e-7);
}


BOOST_AUTO_TEST_CASE ( hubbard_upperTriagonal ) {
    
    Eigen::VectorXd triagonal_test(6);
    triagonal_test << 1, 2, 3, 4, 5, 6;
    auto hubbard_ham_par = GQCP::hubbardTriagonalLattice(triagonal_test);

    Eigen::MatrixXd h_ref (3, 3);
    h_ref << 0, 2, 3,
             2, 0, 5,
             3, 5, 0;

    Eigen::Tensor<double, 4> g_ref (3, 3, 3, 3);
    g_ref.setZero();
    g_ref(0,0,0,0) = 1;
    g_ref(1,1,1,1) = 4;
    g_ref(2,2,2,2) = 6;

    BOOST_CHECK(h_ref.isApprox(hubbard_ham_par.get_h().get_matrix_representation()));
    BOOST_CHECK(cpputil::linalg::areEqual(g_ref, hubbard_ham_par.get_g().get_matrix_representation(), 1.0e-12));
    
    Eigen::VectorXd triagonal_test_faulty(5);
    triagonal_test_faulty << 1, 2, 3, 4, 5;
    
    BOOST_CHECK_THROW(GQCP::hubbardTriagonalLattice(triagonal_test_faulty), std::invalid_argument);
}