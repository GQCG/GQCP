#define BOOST_TEST_MODULE "DOCISolver"


#include "CISolver/CISolver.hpp"
#include "HamiltonianBuilder/DOCI.hpp"
#include "HamiltonianParameters/HamiltonianParameters_constructors.hpp"
#include "RHF/PlainRHFSCFSolver.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain



BOOST_AUTO_TEST_CASE ( DOCI_beh_cation_klaas_dense ) {

    // Klaas' reference DOCI energy for BeH+ (obtained through Caitlin)
    double reference_doci_energy = -14.8782216937;

    // Do a DOCI calculation based on a given FCIDUMP file
    // Create the Hamiltonian Parameters
    auto ham_par = GQCG::readFCIDUMPFile("../tests/data/beh_cation_631g_caitlin.FCIDUMP");

    // The species contains 4 electrons and 16 basis functions, this requires a single Fock Space of 16 orbitals and 2 electrons
    GQCG::FockSpace fock_space (16, 2);  // dim = 120

    // Create the DOCI module
    GQCG::DOCI doci (fock_space);

    // Solve
    GQCG::CISolver ci_solver (doci, ham_par);
    numopt::eigenproblem::DenseSolverOptions solver_options;
    ci_solver.solve(solver_options);

    // Retrieve the eigenvalues
    auto doci_eigenvalue = ci_solver.get_eigenpair().get_eigenvalue();

    // Calculate the total energy
    double internuclear_repulsion_energy = 1.5900757460937498e+00;  // this comes straight out of the FCIDUMP file
    double test_doci_energy = doci_eigenvalue + internuclear_repulsion_energy;

    BOOST_CHECK(std::abs(test_doci_energy - (reference_doci_energy)) < 1.0e-9);
}


BOOST_AUTO_TEST_CASE ( DOCI_lih_klaas_dense ) {

    // Klaas' reference DOCI energy for LiH (obtained through Caitlin)
    double reference_doci_energy = -8.0029560313;

    // Do a DOCI calculation based on a given FCIDUMP file
    // Create the Hamiltonian Parameters
    auto ham_par = GQCG::readFCIDUMPFile("../tests/data/lih_631g_caitlin.FCIDUMP");

    // The species contains 4 electrons and 16 basis functions, this requires a single Fock Space of 16 orbitals and 2 electrons
    GQCG::FockSpace fock_space (16, 2);  // dim = 120

    // Create the DOCI module
    GQCG::DOCI doci (fock_space);

    // Solve
    GQCG::CISolver ci_solver (doci, ham_par);
    numopt::eigenproblem::DenseSolverOptions solver_options;
    ci_solver.solve(solver_options);

    // Retrieve the eigenvalues
    auto doci_eigenvalue = ci_solver.get_eigenpair().get_eigenvalue();

    // Calculate the total energy
    double internuclear_repulsion_energy = 9.6074293445896852e-01;  // this comes straight out of the FCIDUMP file
    double test_doci_energy = doci_eigenvalue + internuclear_repulsion_energy;

    BOOST_CHECK(std::abs(test_doci_energy - (reference_doci_energy)) < 1.0e-9);
}


BOOST_AUTO_TEST_CASE ( DOCI_li2_klaas_dense ) {

    // Klaas' reference DOCI energy for Li2
    double reference_doci_energy = -15.1153976060;

    // Do a DOCI calculation based on a given FCIDUMP file
    // Create the Hamiltonian Parameters
    auto ham_par = GQCG::readFCIDUMPFile("../tests/data/li2_321g_klaas.FCIDUMP");

    // The species contains 4 electrons and 16 basis functions, this requires a single Fock Space of 16 orbitals and 2 electrons
    GQCG::FockSpace fock_space (18, 3);  // dim = 816

    // Create the DOCI module
    GQCG::DOCI doci (fock_space);

    // Solve
    GQCG::CISolver ci_solver (doci, ham_par);
    numopt::eigenproblem::DenseSolverOptions solver_options;
    ci_solver.solve(solver_options);

    // Retrieve the eigenvalues
    auto doci_eigenvalue = ci_solver.get_eigenpair().get_eigenvalue();

    // Calculate the total energy
    double internuclear_repulsion_energy =  3.0036546888874875e+00;  // this comes straight out of the FCIDUMP file
    double test_doci_energy = doci_eigenvalue + internuclear_repulsion_energy;

    BOOST_CHECK(std::abs(test_doci_energy - (reference_doci_energy)) < 1.0e-9);
}


BOOST_AUTO_TEST_CASE ( DOCI_h2o_sto3g_klaas_Davidson ) {

    // Klaas' reference DOCI energy for H2O@STO-3G
    double reference_doci_energy = -74.9671366903;


    // Do a DOCI calculation based on a given FCIDUMP file
    // Create the Hamiltonian Parameters
    auto ham_par = GQCG::readFCIDUMPFile("../tests/data/h2o_sto3g_klaas.FCIDUMP");

    // The species contains 10 electrons and 7 basis functions, this requires a single Fock Space of 7 orbitals and 5 electrons
    GQCG::FockSpace fock_space (7, 5);  // dim = 21

    // Create the DOCI module
    GQCG::DOCI doci (fock_space);

    // Solve
    GQCG::CISolver ci_solver (doci, ham_par);
    Eigen::VectorXd initial_g = fock_space.hartreeFockGuess();
    numopt::eigenproblem::DavidsonSolverOptions solver_options (initial_g);
    ci_solver.solve(solver_options);

    // Retrieve the eigenvalues
    auto doci_eigenvalue = ci_solver.get_eigenpair().get_eigenvalue();

    // Calculate the total energy
    double internuclear_repulsion_energy =  9.7794061444134091E+00;  // this comes straight out of the FCIDUMP file
    double test_doci_energy = doci_eigenvalue + internuclear_repulsion_energy;

    BOOST_CHECK(std::abs(test_doci_energy - (reference_doci_energy)) < 1.0e-9);
}

