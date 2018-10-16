#define BOOST_TEST_MODULE "DenseDOCISolver"


#include "CISolver/CISolver.hpp"
#include "HamiltonianBuilder/DOCI.hpp"
#include "HamiltonianParameters/HamiltonianParameters_constructors.hpp"
#include "RHF/PlainRHFSCFSolver.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain



BOOST_AUTO_TEST_CASE ( test_random_rotation_diagonal_dense_fci ) {

    // Check if a random rotation has no effect on the sum of the diagonal elements

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

    GQCG::FockSpaceProduct fock_space (K, h2o.get_N()/2, h2o.get_N()/2);  // dim = 2

    // Create the FCI module
    GQCG::FCI fci (fock_space);

    Eigen::VectorXd diagonal1 = fci.calculateDiagonal(mol_ham_par);

    // Get a random unitary matrix by diagonalizing a random symmetric matrix
    Eigen::MatrixXd A_random = Eigen::MatrixXd::Random(K, K);
    Eigen::MatrixXd A_symmetric = A_random + A_random.transpose();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> unitary_solver (A_symmetric);
    Eigen::MatrixXd U_random = unitary_solver.eigenvectors();

    // Rotate the hampar using the random unitary matrix
    mol_ham_par.rotate(U_random);

    Eigen::VectorXd diagonal2 = fci.calculateDiagonal(mol_ham_par);

    // Solve Dense
    numopt::eigenproblem::DenseSolverOptions dense_solver_options;
    ci_solver.solve(dense_solver_options);

    // Retrieve the eigenvalues
    auto fci_dense_eigenvalue = ci_solver.get_eigenpair().get_eigenvalue();

    BOOST_CHECK(std::abs(diagonal1.sum() - diagonal2.sum()) < 1.0e-10);
}


BOOST_AUTO_TEST_CASE ( FCI_H2_Cristina_dense ) {

    // Cristina's H2 FCI energy/OO-DOCI energy
    double reference_fci_energy = -1.1651486697;

    // Create a Molecule and an AOBasis
    GQCG::Molecule h2 ("../tests/data/h2_cristina.xyz");
    auto ao_basis = std::make_shared<GQCG::AOBasis>(h2, "6-31g**");

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

    // Solve Dense
    numopt::eigenproblem::DenseSolverOptions dense_solver_options;
    ci_solver.solve(dense_solver_options);

    // Retrieve the eigenvalues
    auto test_fci_energy = ci_solver.get_eigenpair().get_eigenvalue();

    // Calculate the total FCI energy
    double internuclear_repulsion_energy = h2.calculateInternuclearRepulsionEnergy();
    double test_fci_energy = fci.get_eigenvalue() + internuclear_repulsion_energy;

    BOOST_CHECK(std::abs(test_fci_energy - (reference_fci_energy)) < 1.0e-06);
}


BOOST_AUTO_TEST_CASE ( FCI_H2O_Psi4_GAMESS_dense ) {

    // Psi4 and GAMESS' FCI energy
    double reference_fci_energy = -75.0129803939602;

    // Create a Molecule and an AOBasis
    GQCG::Molecule h2o ("../tests/data/h2o_Psi4_GAMESS.xyz");
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

    // Solve Dense
    numopt::eigenproblem::DenseSolverOptions dense_solver_options;
    ci_solver.solve(dense_solver_options);

    // Retrieve the eigenvalues
    auto test_fci_energy = ci_solver.get_eigenpair().get_eigenvalue();

    // Calculate the total FCI energy
    double internuclear_repulsion_energy = h2o.calculateInternuclearRepulsionEnergy();
    double test_fci_energy = fci.get_eigenvalue() + internuclear_repulsion_energy;

    BOOST_CHECK(std::abs(test_fci_energy - (reference_fci_energy)) < 1.0e-06);
}


BOOST_AUTO_TEST_CASE ( FCI_He_Cristina_dense ) {

    // Cristina's He FCI energy
    double reference_fci_energy = -2.902533599;

    // Create a Molecule and an AOBasis
    GQCG::Molecule he ("../tests/data/h2o_Psi4_GAMESS.xyz");
    auto ao_basis = std::make_shared<GQCG::AOBasis>(he, "aug-cc-pVQZ");

    // Create the molecular Hamiltonian parameters for this molecule and basis
    auto mol_ham_par = GQCG::constructMolecularHamiltonianParameters(ao_basis);
    auto K = mol_ham_par.get_h().get_dim();

    // Create a plain RHF SCF solver and solve the SCF equations
    GQCG::PlainRHFSCFSolver plain_scf_solver (mol_ham_par, he);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    // Transform the ham_par
    mol_ham_par.transform(rhf.get_C());

    GQCG::FockSpaceProduct fock_space (K, he.get_N()/2, he.get_N()/2);  // dim = 2116

    // Create the FCI module
    GQCG::FCI fci (fock_space);
    GQCG::CISolver ci_solver (fci, mol_ham_par);

    // Solve Dense
    numopt::eigenproblem::DenseSolverOptions dense_solver_options;
    ci_solver.solve(dense_solver_options);

    // Retrieve the eigenvalues
    auto test_fci_energy = ci_solver.get_eigenpair().get_eigenvalue();

    // Calculate the total FCI energy
    double internuclear_repulsion_energy = he.calculateInternuclearRepulsionEnergy();
    double test_fci_energy = fci.get_eigenvalue() + internuclear_repulsion_energy;

    BOOST_CHECK(std::abs(test_fci_energy - (reference_fci_energy)) < 1.0e-06);
}