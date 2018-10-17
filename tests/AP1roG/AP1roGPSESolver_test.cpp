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
    size_t N = 2;  // number of electrons for H2
    size_t N_P = N/2;  // number of electron pairs for H2
    auto ao_basis = std::make_shared<GQCG::AOBasis>(h2, "STO-3G");
    auto mol_ham_par = GQCG::constructMolecularHamiltonianParameters(ao_basis);
    GQCG::AP1roGPSESolver ap1rog_pse_solver (N_P, mol_ham_par);
}


BOOST_AUTO_TEST_CASE ( constructor_molecule ) {

    // Test a correct constructor
    // Check if we can also pass a molecule object to the constructor
    GQCG::Molecule h2 ("../tests/data/h2_szabo.xyz");
    auto ao_basis = std::make_shared<GQCG::AOBasis>(h2, "STO-3G");
    auto mol_ham_par = GQCG::constructMolecularHamiltonianParameters(ao_basis);
    GQCG::AP1roGPSESolver ap1rog_pse_solver (h2, mol_ham_par);

    // Test a faulty constructor
    GQCG::Molecule h2_cation ("../tests/data/h2_szabo.xyz", +1);
    BOOST_CHECK_THROW(GQCG::AP1roGPSESolver(h2_cation, mol_ham_par), std::invalid_argument);  // we can use the same Hamiltonian parameters for molecule and ion
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
