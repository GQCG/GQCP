#define BOOST_TEST_MODULE "AP1roG"


#include "AP1roG/AP1roGPSESolver.hpp"
#include "HamiltonianParameters/HamiltonianParameters_constructors.hpp"
#include "Molecule.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


BOOST_AUTO_TEST_CASE ( constructor ) {

    // Test a correct constructor
    GQCG::Molecule h2 ("../tests/data/h2_szabo.xyz");
    auto ao_basis = std::make_shared<GQCG::AOBasis>(h2, "STO-3G");
    auto mol_ham_par = GQCG::constructMolecularHamiltonianParameters(ao_basis);
    GQCG::AP1roGPSESolver ap1rog_pse_solver (h2, mol_ham_par);


    // Test a faulty constructor
    GQCG::Molecule h2_ion ("../tests/data/h2_szabo.xyz", +1);
    BOOST_CHECK_THROW(GQCG::AP1roGPSESolver(h2_ion, mol_ham_par), std::invalid_argument);  // we can use the same Hamiltonian parameters for molecule and ion
}


BOOST_AUTO_TEST_CASE ( vector_index ) {

    // Create a real example for testing the vector index
    // In this case we have LiH//6-31G
    // N_P=2, K=6
    GQCG::Molecule lih ("../tests/data/lih_olsens.xyz");
    auto ao_basis = std::make_shared<GQCG::AOBasis>(lih, "6-31G");
    auto mol_ham_par = GQCG::constructMolecularHamiltonianParameters(ao_basis);
    GQCG::AP1roGPSESolver ap1rog_pse_solver (lih, mol_ham_par);


    BOOST_CHECK_EQUAL(ap1rog_pse_solver.vectorIndex(0, 2), 0);
    BOOST_CHECK_EQUAL(ap1rog_pse_solver.vectorIndex(0, 3), 1);
    BOOST_CHECK_EQUAL(ap1rog_pse_solver.vectorIndex(1, 2), 9);
    BOOST_CHECK_EQUAL(ap1rog_pse_solver.vectorIndex(1, 3), 10);

    // Require a throw if i > N_P
    BOOST_REQUIRE_THROW(ap1rog_pse_solver.vectorIndex(3, 3), std::invalid_argument);

    // Require a throw if a < N_P
    BOOST_REQUIRE_THROW(ap1rog_pse_solver.vectorIndex(0, 1), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( matrix_index ) {

    // Create a real example for testing the vector index
    // In this case we have LiH//6-31G
    // N_P=2, K=6
    GQCG::Molecule lih ("../tests/data/lih_olsens.xyz");
    auto ao_basis = std::make_shared<GQCG::AOBasis>(lih, "6-31G");
    auto mol_ham_par = GQCG::constructMolecularHamiltonianParameters(ao_basis);
    GQCG::AP1roGPSESolver ap1rog_pse_solver (lih, mol_ham_par);


    BOOST_CHECK_EQUAL(ap1rog_pse_solver.matrixIndexMajor(0), 0);
    BOOST_CHECK_EQUAL(ap1rog_pse_solver.matrixIndexMajor(1), 0);
    BOOST_CHECK_EQUAL(ap1rog_pse_solver.matrixIndexMajor(9), 1);
    BOOST_CHECK_EQUAL(ap1rog_pse_solver.matrixIndexMajor(10), 1);

    BOOST_CHECK_EQUAL(ap1rog_pse_solver.matrixIndexMinor(0), 2);
    BOOST_CHECK_EQUAL(ap1rog_pse_solver.matrixIndexMinor(1), 3);
    BOOST_CHECK_EQUAL(ap1rog_pse_solver.matrixIndexMinor(4), 6);
    BOOST_CHECK_EQUAL(ap1rog_pse_solver.matrixIndexMinor(5), 7);
}
