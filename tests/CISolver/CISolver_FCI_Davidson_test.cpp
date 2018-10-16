#define BOOST_TEST_MODULE "DavidsonFCISolver"


#include "CISolver/CISolver.hpp"
#include "HamiltonianBuilder/FCI.hpp"
#include "HamiltonianParameters/HamiltonianParameters_constructors.hpp"
#include "RHF/PlainRHFSCFSolver.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain




BOOST_AUTO_TEST_CASE ( FCI_h2_sto3g_dense_vs_Davidson ) {

    // Check if the dense FCI energy is equal to the Davidson (with matvec) FCI energy

    // Create a Molecule and an AOBasis
    GQCG::Molecule h2 ("../tests/data/h2_cristina.xyz");
    auto ao_basis = std::make_shared<GQCG::AOBasis>(h2, "STO-3G");

    // Create the molecular Hamiltonian parameters for this molecule and basis
    auto mol_ham_par = GQCG::constructMolecularHamiltonianParameters(ao_basis);
    auto K = mol_ham_par.get_h().get_dim();

    // Create a plain RHF SCF solver and solve the SCF equations
    GQCG::PlainRHFSCFSolver plain_scf_solver (mol_ham_par, h2);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    // Transform the ham_par
    mol_ham_par.transform(rhf.get_C());

    GQCG::FockSpaceProduct fock_space (K, h2.get_N()/2, h2.get_N()/2);  // dim = 2

    // Create the FCI module
    GQCG::FCI fci (fock_space);
    GQCG::CISolver ci_solver (fci, mol_ham_par);

    // Solve Davidson
    Eigen::VectorXd initial_g = fock_space.HartreeFockExpansion();
    numopt::eigenproblem::DavidsonSolverOptions davidson_solver_options (initial_g);
    ci_solver.solve(davidson_solver_options);

    // Retrieve the eigenvalues
    auto fci_davidson_eigenvalue = ci_solver.get_eigenpair().get_eigenvalue();

    // Solve Dense
    numopt::eigenproblem::DenseSolverOptions dense_solver_options;
    ci_solver.solve(dense_solver_options);

    // Retrieve the eigenvalues
    auto fci_dense_eigenvalue = ci_solver.get_eigenpair().get_eigenvalue();

    BOOST_CHECK(std::abs(fci_dense_eigenvalue - fci_davidson_eigenvalue) < 1.0e-08);
}


BOOST_AUTO_TEST_CASE ( FCI_H2_6_31Gxx_dense_vs_Davidson ) {

    // Check if the dense FCI energy is equal to the Davidson (with matvec) FCI energy

    // Create a Molecule and an AOBasis
    GQCG::Molecule h2 ("../tests/data/h2_cristina.xyz");
    auto ao_basis = std::make_shared<GQCG::AOBasis>(h2, "6-31G**");

    // Create the molecular Hamiltonian parameters for this molecule and basis
    auto mol_ham_par = GQCG::constructMolecularHamiltonianParameters(ao_basis);
    auto K = mol_ham_par.get_h().get_dim();

    // Create a plain RHF SCF solver and solve the SCF equations
    GQCG::PlainRHFSCFSolver plain_scf_solver (mol_ham_par, h2);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    // Transform the ham_par
    mol_ham_par.transform(rhf.get_C());

    GQCG::FockSpaceProduct fock_space (K, h2.get_N()/2, h2.get_N()/2);  // dim = 100

    // Create the FCI module
    GQCG::FCI fci (fock_space);
    GQCG::CISolver ci_solver (fci, mol_ham_par);

    // Solve Davidson
    Eigen::VectorXd initial_g = fock_space.HartreeFockExpansion();
    numopt::eigenproblem::DavidsonSolverOptions davidson_solver_options (initial_g);
    ci_solver.solve(davidson_solver_options);

    // Retrieve the eigenvalues
    auto fci_davidson_eigenvalue = ci_solver.get_eigenpair().get_eigenvalue();

    // Solve Dense
    numopt::eigenproblem::DenseSolverOptions dense_solver_options;
    ci_solver.solve(dense_solver_options);

    // Retrieve the eigenvalues
    auto fci_dense_eigenvalue = ci_solver.get_eigenpair().get_eigenvalue();

    BOOST_CHECK(std::abs(fci_dense_eigenvalue - fci_davidson_eigenvalue) < 1.0e-08);
}


BOOST_AUTO_TEST_CASE ( FCI_H2O_STO_3G_dense_vs_Davidson ) {

    // Check if the dense FCI energy is equal to the Davidson (with matvec) FCI energy

    // Create a Molecule and an AOBasis
    GQCG::Molecule h2o ("../tests/data/h2o.xyz");
    auto ao_basis = std::make_shared<GQCG::AOBasis>(h2o, "STO-3G");

    // Create the molecular Hamiltonian parameters for this molecule and basis
    auto mol_ham_par = GQCG::constructMolecularHamiltonianParameters(ao_basis);
    auto K = mol_ham_par.get_h().get_dim();

    // Create a plain RHF SCF solver and solve the SCF equations
    GQCG::PlainRHFSCFSolver plain_scf_solver (mol_ham_par, h2o);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    // Transform the ham_par
    mol_ham_par.transform(rhf.get_C());

    GQCG::FockSpaceProduct fock_space (K, h2o.get_N()/2, h2o.get_N()/2);  // dim = 441

    // Create the FCI module
    GQCG::FCI fci (fock_space);
    GQCG::CISolver ci_solver (fci, mol_ham_par);

    // Solve Davidson
    Eigen::VectorXd initial_g = fock_space.HartreeFockExpansion();
    numopt::eigenproblem::DavidsonSolverOptions davidson_solver_options (initial_g);
    ci_solver.solve(davidson_solver_options);

    // Retrieve the eigenvalues
    auto fci_davidson_eigenvalue = ci_solver.get_eigenpair().get_eigenvalue();

    // Solve Dense
    numopt::eigenproblem::DenseSolverOptions dense_solver_options;
    ci_solver.solve(dense_solver_options);

    // Retrieve the eigenvalues
    auto fci_dense_eigenvalue = ci_solver.get_eigenpair().get_eigenvalue();

    BOOST_CHECK(std::abs(fci_dense_eigenvalue - fci_davidson_eigenvalue) < 1.0e-08);
}