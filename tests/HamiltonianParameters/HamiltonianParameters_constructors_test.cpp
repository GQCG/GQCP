#define BOOST_TEST_MODULE "HamiltonianParameters_constructors"


#include "HamiltonianParameters/HamiltonianParameters_constructors.hpp"


#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain



BOOST_AUTO_TEST_CASE ( constructMolecularHamiltonianParameters ) {

    // Set up a basis
    GQCG::Molecule water ("../tests/data/h2o.xyz");
    auto ao_basis_sptr = std::make_shared<GQCG::AOBasis>(water, "STO-3G");

    
    // Check if we can construct the molecular Hamiltonian parameters
    auto mol_ham_par = GQCG::constructMolecularHamiltonianParameters(ao_basis_sptr);
}


BOOST_AUTO_TEST_CASE ( FCIDUMP_reader ) {
    
    auto fcidump_ham_par = GQCG::readFCIDUMPFile("../tests/data/beh_cation_631g_caitlin.FCIDUMP");
    
    // Check if the one-electron integrals are read in correctly from a previous implementation
    Eigen::MatrixXd h_SO = fcidump_ham_par.get_h().get_matrix_representation();
    
    BOOST_CHECK(std::abs(h_SO(0,0) - (-8.34082)) < 1.0e-5);
    BOOST_CHECK(std::abs(h_SO(5,1) - 0.381418) < 1.0e-6);
    BOOST_CHECK(std::abs(h_SO(14,0) - 0.163205) < 1.0e-6);
    BOOST_CHECK(std::abs(h_SO(13,6) - (-5.53204e-16)) < 1.0e-16);
    BOOST_CHECK(std::abs(h_SO(15,11) - (-0.110721)) < 1.0e-6);
    
    
    // Check if the two-electron integrals are read in correctly from a previous implementation
    Eigen::Tensor<double, 4> g_SO = fcidump_ham_par.get_g().get_matrix_representation();
    
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
    auto fcidump_ham_par = GQCG::readFCIDUMPFile("../tests/data/h2_psi4_horton.FCIDUMP");
    
    Eigen::Tensor<double, 4> g_SO = fcidump_ham_par.get_g().get_matrix_representation();
    BOOST_CHECK(std::abs(g_SO(6,5,1,0) - 0.0533584656) <  1.0e-7);
}
