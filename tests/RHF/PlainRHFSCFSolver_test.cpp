#define BOOST_TEST_MODULE "PlainRHFSCFSolver"

#include "RHF/PlainRHFSCFSolver.hpp"
#include "HamiltonianParameters/HamiltonianParameters_constructors.hpp"

#include <cpputil.hpp>

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain



BOOST_AUTO_TEST_CASE ( h2_sto3g_szabo_plain ) {

    // In this test case, we will follow section 3.5.2 in Szabo.
    double ref_total_energy = -1.1167;


    // Create a Molecule and an AOBasis
    GQCP::Molecule h2 ("../tests/data/h2_szabo.xyz");
    auto ao_basis = std::make_shared<GQCP::AOBasis>(h2, "STO-3G");

    // Create the molecular Hamiltonian parameters for this molecule and basis
    auto mol_ham_par = GQCP::constructMolecularHamiltonianParameters(ao_basis);


    // Create a plain RHF SCF solver and solve the SCF equations
    GQCP::PlainRHFSCFSolver plain_scf_solver (mol_ham_par, h2);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();


    // Check the total energy
    double total_energy = rhf.get_electronic_energy() + h2.calculateInternuclearRepulsionEnergy();
    BOOST_CHECK(std::abs(total_energy - ref_total_energy) < 1.0e-04);
}


BOOST_AUTO_TEST_CASE ( h2o_sto3g_horton_plain ) {

    // We have some reference data from horton
    double ref_total_energy = -74.942080055631;

    Eigen::VectorXd ref_orbital_energies (7);  // the STO-3G basisset has 7 basis functions for water
    ref_orbital_energies << -20.26289322, -1.20969863, -0.54796582, -0.43652631, -0.38758791, 0.47762043, 0.5881361;

    Eigen::MatrixXd ref_C (7, 7);
    ref_C << -9.94434594e-01, -2.39158997e-01,  3.61117086e-17, -9.36837259e-02,  3.73303682e-31, -1.11639152e-01, -9.04958229e-17,
             -2.40970260e-02,  8.85736467e-01, -1.62817254e-16,  4.79589270e-01, -1.93821120e-30,  6.69575233e-01,  5.16088339e-16,
              1.59542752e-18,  5.29309704e-17, -6.07288675e-01, -1.49717339e-16,  8.94470461e-17, -8.85143477e-16,  9.19231270e-01,
             -3.16155527e-03,  8.58957413e-02,  2.89059171e-16, -7.47426286e-01,  2.81871324e-30,  7.38494291e-01,  6.90314422e-16,
              6.65079968e-35,  1.16150362e-32, -2.22044605e-16, -4.06685146e-30, -1.00000000e+00, -1.78495825e-31,  2.22044605e-16,
              4.59373756e-03,  1.44038811e-01, -4.52995183e-01, -3.29475784e-01,  2.16823939e-16, -7.09847234e-01, -7.32462496e-01,
              4.59373756e-03,  1.44038811e-01,  4.52995183e-01, -3.29475784e-01, -2.16823939e-16, -7.09847234e-01,  7.32462496e-01;


    // Do our own RHF calculation
    GQCP::Molecule water ("../tests/data/h2o.xyz");
    auto ao_basis = std::make_shared<GQCP::AOBasis>(water, "STO-3G");
    auto mol_ham_par = GQCP::constructMolecularHamiltonianParameters(ao_basis);

    GQCP::PlainRHFSCFSolver plain_scf_solver (mol_ham_par, water);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    // Check the total energy
    double total_energy = rhf.get_electronic_energy() + water.calculateInternuclearRepulsionEnergy();
    BOOST_CHECK(std::abs(total_energy - ref_total_energy) < 1.0e-04);

    // Check the calculated results with the reference
    BOOST_CHECK(std::abs(total_energy - ref_total_energy) < 1.0e-06);
    BOOST_CHECK(cpputil::linalg::areEqualEigenvalues(ref_orbital_energies, rhf.get_orbital_energies(), 1.0e-06));
    BOOST_CHECK(cpputil::linalg::areEqualSetsOfEigenvectors(ref_C, rhf.get_C(), 1.0e-05));
}


BOOST_AUTO_TEST_CASE ( crawdad_h2o_sto3g_plain ) {

    // This example is taken from (http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project3), but the input .xyz-file was converted to Angstrom.
    double ref_total_energy = -74.9420799281920;


    // Do our own RHF calculation
    GQCP::Molecule water ("../tests/data/h2o_crawdad.xyz");
    // Check if the internuclear distance between O and H is really 1.1 A (= 2.07869 bohr), as specified in the text
    BOOST_REQUIRE(std::abs(water.calculateInternuclearDistance(0, 1) - 2.07869) < 1.0e-4);

    auto ao_basis = std::make_shared<GQCP::AOBasis>(water, "STO-3G");
    auto mol_ham_par = GQCP::constructMolecularHamiltonianParameters(ao_basis);

    GQCP::PlainRHFSCFSolver plain_scf_solver (mol_ham_par, water);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();


    // Check the total energy
    double total_energy = rhf.get_electronic_energy() + water.calculateInternuclearRepulsionEnergy();
    BOOST_CHECK(std::abs(total_energy - ref_total_energy) < 1.0e-06);
}


BOOST_AUTO_TEST_CASE ( crawdad_ch4_sto3g_plain ) {

    // This example is taken from (http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project3), but the input .xyz-file was converted to Angstrom.
    double ref_total_energy = -39.726850324347;

    // Do our own RHF calculation
    GQCP::Molecule methane ("../tests/data/ch4_crawdad.xyz");
    // Check if the internuclear distance between C and H is really around 2.05 bohr, which is the bond distance Wikipedia (108.7 pm) specifies
    BOOST_CHECK(std::abs(methane.calculateInternuclearDistance(0, 1) - 2.05) < 1.0e-1);

    auto ao_basis = std::make_shared<GQCP::AOBasis>(methane, "STO-3G");
    auto mol_ham_par = GQCP::constructMolecularHamiltonianParameters(ao_basis);

    GQCP::PlainRHFSCFSolver plain_scf_solver (mol_ham_par, methane);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    // Check the total energy
    double total_energy = rhf.get_electronic_energy() + methane.calculateInternuclearRepulsionEnergy();
    BOOST_CHECK(std::abs(total_energy - ref_total_energy) < 1.0e-06);
}


BOOST_AUTO_TEST_CASE ( h2_631gdp_plain ) {

    // We have some reference data from olsens: H2@RHF//6-31G** orbitals
    double ref_electronic_energy = -1.84444667247;


    // Do our own RHF calculation
    GQCP::Molecule h2 ("../tests/data/h2_olsens.xyz");
    auto ao_basis = std::make_shared<GQCP::AOBasis>(h2, "6-31g**");
    auto mol_ham_par = GQCP::constructMolecularHamiltonianParameters(ao_basis);

    GQCP::PlainRHFSCFSolver plain_scf_solver (mol_ham_par, h2);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    // Check the electronic energy
    BOOST_CHECK(std::abs(rhf.get_electronic_energy() - ref_electronic_energy) < 1.0e-06);
}
