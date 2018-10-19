#define BOOST_TEST_MODULE "RMP2"

#include "RMP2.hpp"

#include "HamiltonianParameters/HamiltonianParameters_constructors.hpp"
#include "RHF/PlainRHFSCFSolver.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain



BOOST_AUTO_TEST_CASE ( crawdad_sto3g_water ) {

    // Get the reference data from crawdad (http://sirius.chem.vt.edu/~crawdad/programming/project4/h2o_sto3g/output.txt)
    double ref_energy_correction = -0.049149636120;


    // Create molecular Hamiltonian parameters in the RHF basis
    GQCP::Molecule water ("../tests/data/h2o_crawdad.xyz");
    auto ao_basis = std::make_shared<GQCP::AOBasis>(water, "STO-3G");
    auto ao_mol_ham_par = GQCP::constructMolecularHamiltonianParameters(ao_basis);

    GQCP::PlainRHFSCFSolver plain_scf_solver (ao_mol_ham_par, water);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();
    auto mol_ham_par = GQCP::HamiltonianParameters(ao_mol_ham_par, rhf.get_C());


    // Check if the RMP2 correction is correct
    double energy_correction = GQCP::calculateRMP2EnergyCorrection(mol_ham_par, water, rhf);
    BOOST_CHECK(std::abs(energy_correction - ref_energy_correction) < 1.0e-08);
}


BOOST_AUTO_TEST_CASE ( crawdad_sto3g_methane ) {

    // Get the reference data from crawdad (http://sirius.chem.vt.edu/~crawdad/programming/project4/ch4_sto3g/output.txt)
    double ref_energy_correction = -0.056046676165;

    // Create molecular Hamiltonian parameters in the RHF basis
    GQCP::Molecule methane ("../tests/data/ch4_crawdad.xyz");
    auto ao_basis = std::make_shared<GQCP::AOBasis>(methane, "STO-3G");
    auto ao_mol_ham_par = GQCP::constructMolecularHamiltonianParameters(ao_basis);

    GQCP::PlainRHFSCFSolver plain_scf_solver (ao_mol_ham_par, methane);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();
    auto mol_ham_par = GQCP::HamiltonianParameters(ao_mol_ham_par, rhf.get_C());


    // Check if the RMP2 correction is correct
    double energy_correction = GQCP::calculateRMP2EnergyCorrection(mol_ham_par, methane, rhf);
    BOOST_CHECK(std::abs(energy_correction - ref_energy_correction) < 1.0e-08);
}

