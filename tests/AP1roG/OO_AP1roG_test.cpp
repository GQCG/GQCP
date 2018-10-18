#define BOOST_TEST_MODULE "OO-AP1roG"

#include "AP1roG/AP1roGJacobiOrbitalOptimizer.hpp"
#include "AP1roG/AP1roGPSESolver.hpp"
#include "HamiltonianParameters/HamiltonianParameters_constructors.hpp"
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
    GQCG::AP1roGJacobiOrbitalOptimizer ap1rog_orbital_optimizer (N_P, mol_ham_par);
}


BOOST_AUTO_TEST_CASE ( constructor_molecule ) {

    // Test a correct constructor
    // Check if we can also pass a molecule object to the constructor
    GQCG::Molecule h2 ("../tests/data/h2_szabo.xyz");
    auto ao_basis = std::make_shared<GQCG::AOBasis>(h2, "STO-3G");
    auto mol_ham_par = GQCG::constructMolecularHamiltonianParameters(ao_basis);
    GQCG::AP1roGJacobiOrbitalOptimizer ap1rog_orbital_optimizer (h2, mol_ham_par);


    // Test a faulty constructor
    GQCG::Molecule h2_cation ("../tests/data/h2_szabo.xyz", +1);
    BOOST_CHECK_THROW(GQCG::AP1roGJacobiOrbitalOptimizer(h2_cation, mol_ham_par), std::invalid_argument);  // we can use the same Hamiltonian parameters for molecule and ion
}


BOOST_AUTO_TEST_CASE ( lih_6_31G_calculateEnergyAfterRotation ) {

    // We have implemented a formula to calculate the rotated AP1roG energy directly, but we have to test it
    // It should be equal to the energy we obtain by rotating the one- and two-electron integrals first


    // Construct the molecular Hamiltonian parameters in the RHF basis
    GQCG::Molecule lih ("../tests/data/lih_olsens.xyz");
    auto ao_basis = std::make_shared<GQCG::AOBasis>(lih, "6-31G");
    auto ao_mol_ham_par = GQCG::constructMolecularHamiltonianParameters(ao_basis);

    GQCG::PlainRHFSCFSolver plain_scf_solver (ao_mol_ham_par, lih);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();
    auto mol_ham_par = GQCG::HamiltonianParameters(ao_mol_ham_par, rhf.get_C());


    // Loop over all possible Jacobi pairs for a given (random) angle and check if the analytical result matches the numerical result
    double theta = 56.71;
    size_t K = ao_basis.get()->get_number_of_basis_functions();

    for (size_t q = 0; q < K; q++) {  // p and q loop over spatial orbitals
        for (size_t p = q + 1; p < K; p++) {  // p > q
            GQCG::JacobiRotationParameters jacobi_rot_par {p, q, theta};

            // Construct a new AP1roG instance, since AP1roG.calculateEnergyAfterRotation overwrites this->so_basis
            // AP1roG.calculateEnergyAfterRotation is a function that is only used in testing
            GQCG::AP1roGPSESolver pse_solver (lih, mol_ham_par);
            pse_solver.solve();
            auto G = pse_solver.get_solution().get_geminal_coefficients();


            // Calculate the analytical energy after rotation
            GQCG::AP1roGJacobiOrbitalOptimizer orbital_optimizer (lih, mol_ham_par);
            orbital_optimizer.calculateJacobiCoefficients(p, q, G);
            double E_rotated_analytical = orbital_optimizer.calculateEnergyAfterJacobiRotation(jacobi_rot_par, G);


            // Calculate the energy after a numerical rotation (using a Jacobi rotation matrix)
            auto mol_ham_par_copy = mol_ham_par;
            mol_ham_par_copy.rotate(jacobi_rot_par);
            double E_rotated_numerical = GQCG::calculateAP1roGEnergy(G, mol_ham_par_copy);


            BOOST_CHECK(std::abs(E_rotated_analytical - E_rotated_numerical) < 1.0e-08);
        }
    }
}


BOOST_AUTO_TEST_CASE ( lih_6_31G_orbitalOptimize ) {

    // Construct the molecular Hamiltonian parameters in the RHF basis
    GQCG::Molecule lih ("../tests/data/lih_olsens.xyz");
    auto ao_basis = std::make_shared<GQCG::AOBasis>(lih, "6-31G");
    auto ao_mol_ham_par = GQCG::constructMolecularHamiltonianParameters(ao_basis);

    GQCG::PlainRHFSCFSolver plain_scf_solver (ao_mol_ham_par, lih);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();
    auto mol_ham_par = GQCG::HamiltonianParameters(ao_mol_ham_par, rhf.get_C());


    // Get the initial AP1roG energy
    GQCG::AP1roGPSESolver pse_solver (lih, mol_ham_par);
    pse_solver.solve();
    double initial_energy = pse_solver.get_solution().get_electronic_energy();


    // Do an AP1roG orbital optimization using Jacobi rotations
    GQCG::AP1roGJacobiOrbitalOptimizer orbital_optimizer (lih, mol_ham_par);
    orbital_optimizer.solve();
    double optimized_energy = orbital_optimizer.get_solution().get_electronic_energy();


    // We don't have reference data, so all we can do is check if orbital optimization lowers the energy
    BOOST_CHECK(optimized_energy < initial_energy);
}
