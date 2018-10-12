#define BOOST_TEST_MODULE "AP1roGPSESolver"


#include "AP1roG/AP1roGPSESolver.hpp"

#include "HamiltonianParameters/HamiltonianParameters_constructors.hpp"
#include "Molecule.hpp"
#include "RHF/PlainRHFSCFSolver.hpp"

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


BOOST_AUTO_TEST_CASE ( h2_631gdp ) {

    // We have some reference data from olsens: H2 with HF/6-31G** orbitals
    //  AP1roG energy: -1.8696828608304892
    //  AP1roG coefficients: [-0.05949796, -0.05454253, -0.03709503, -0.02899231, -0.02899231, -0.01317386, -0.00852702, -0.00852702, -0.00517996]

    double ref_ap1rog_energy = -1.8696828608304892;
    Eigen::VectorXd ref_ap1rog_coefficients (9);
    ref_ap1rog_coefficients << -0.05949796, -0.05454253, -0.03709503, -0.02899231, -0.02899231, -0.01317386, -0.00852702, -0.00852702, -0.00517996;


    // 1. Do an RHF calculation
    GQCG::Molecule h2 ("../tests/data/h2_olsens.xyz");
    auto ao_basis = std::make_shared<GQCG::AOBasis>(h2, "6-31G**");
    auto ao_mol_ham_par = GQCG::constructMolecularHamiltonianParameters(ao_basis);

    GQCG::PlainRHFSCFSolver plain_scf_solver (ao_mol_ham_par, h2);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();



    // 2. Solve the AP1roG pSE equations
    auto mol_ham_par = GQCG::HamiltonianParameters(ao_mol_ham_par, rhf.get_C());
    GQCG::AP1roGPSESolver ap1rog_pse_solver (h2, mol_ham_par);
    ap1rog_pse_solver.solve();

    auto ap1rog = ap1rog_pse_solver.get_solution();
    double electronic_energy = ap1rog.get_electronic_energy();
    Eigen::VectorXd ap1rog_coefficients = ap1rog.get_geminal_coefficients().asVector();


    // Check the energy and coefficients
    BOOST_CHECK(std::abs(electronic_energy - ref_ap1rog_energy) < 1.0e-05);

    for (size_t i = 0; i < 9; i++) {
        BOOST_CHECK(std::abs(ap1rog_coefficients(i) - ref_ap1rog_coefficients(i)) < 1.0e-05);
    }
}
