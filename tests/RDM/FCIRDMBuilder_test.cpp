#define BOOST_TEST_MODULE "FCI_RDM_test"



#include "RDM/FCIRDMBuilder.hpp"

#include "CISolver/CISolver.hpp"
#include "HamiltonianBuilder/FCI.hpp"
#include "HamiltonianParameters/HamiltonianParameters_constructors.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>


BOOST_AUTO_TEST_CASE ( lih_1RDM_trace ) {

    // Test if the trace of the 1-RDM gives N

    // Get the 1-RDM from FCI
    size_t N_a = 5;
    size_t N_b = 5;
    
    auto ham_par = GQCG::readFCIDUMPFile("../tests/data/h2o_Psi4_GAMESS.xyz");
    size_t K = ham_par.get_K();  // SO 7

    GQCG::FockSpaceProduct fock_space (K, N_a, N_b);  // dim = 441
    GQCG::FCI fci (fock_space);

    // Specify solver options and solve the eigenvalue problem
    // Solve the dense FCI eigenvalue problem
    GQCG::CISolver ci_solver (fci, ham_par);
    numopt::eigenproblem::DenseSolverOptions solver_options;
    ci_solver.solve(solver_options);

    Eigen::VectorXd coef = ci_solver.get_eigenpair().get_eigenvector();

    // Check if the FCI 1-RDM has the proper trace.
    GQCG::FCIRDMBuilder fci_rdm (fock_space);
    GQCG::OneRDMs one_rdms = fci_rdm.calculate1RDMs(coef);

    BOOST_CHECK(std::abs(one_rdms.one_rdm.trace() - N) < 1.0e-12);
}


BOOST_AUTO_TEST_CASE ( lih_2RDM_trace ) {

    // Test if the trace of the 2-RDM (d_ppqq) gives N(N-1)


    // Get the 2-RDM from FCI
    size_t N = 4;  // 4 electrons
    size_t K = 16;  // 16 SOs
    auto ham_par = GQCG::readFCIDUMPFile("../tests/data/lih_631g_caitlin.FCIDUMP");

    GQCG::FockSpaceProduct fock_space (16, 2);  // dim = 120
    GQCG::FCI fci (fock_space);

    // Specify solver options and solve the eigenvalue problem
    // Solve the dense FCI eigenvalue problem
    GQCG::CISolver ci_solver (fci, ham_par);
    numopt::eigenproblem::DenseSolverOptions solver_options;
    ci_solver.solve(solver_options);

    Eigen::VectorXd coef = ci_solver.get_eigenpair().get_eigenvector();

    // Check if the 2-RDM has the proper trace.
    GQCG::FCIRDMBuilder fci_rdm (fock_space);
    GQCG::TwoRDMs two_rdms = fci_rdm.calculate2RDMs(coef);

    BOOST_CHECK(std::abs(two_rdms.two_rdm.trace() - N*(N-1)) < 1.0e-12);
}


BOOST_AUTO_TEST_CASE ( lih_1RDM_2RDM_trace_FCI ) {

    // Test if the relevant 2-RDM trace gives the 1-RDM for FCI


    // Get the 1- and 2-RDMs from FCI
    size_t N = 4;  // 4 electrons
    size_t K = 16;  // 16 SOs
    auto ham_par = GQCG::readFCIDUMPFile("../tests/data/lih_631g_caitlin.FCIDUMP");

    GQCG::FockSpaceProduct fock_space (16, 2);  // dim = 120
    GQCG::FCI fci (fock_space);

    // Specify solver options and solve the eigenvalue problem
    // Solve the dense FCI eigenvalue problem
    GQCG::CISolver ci_solver (fci, ham_par);
    numopt::eigenproblem::DenseSolverOptions solver_options;
    ci_solver.solve(solver_options);

    Eigen::VectorXd coef = ci_solver.get_eigenpair().get_eigenvector();

    // Check if the 2-RDM contraction matches the reduction.
    GQCG::FCIRDMBuilder fci_rdm (fock_space);
    GQCG::TwoRDMs two_rdms = fci_rdm.calculate2RDMs(coef);
    GQCG::OneRDMs one_rdms = fci_rdm.calculate1RDMs(coef);


    Eigen::MatrixXd D_from_reduction = (1.0/(N-1)) * two_rdms.two_rdm.reduce();
    BOOST_CHECK(one_rdms.one_rdm.get_matrix_representation().isApprox(D_from_reduction, 1.0e-12));
}


BOOST_AUTO_TEST_CASE ( lih_energy_RDM_contraction_FCI ) {

    // Test if the contraction of the 1- and 2-RDMs with the one- and two-electron integrals gives the FCI energy

    size_t N = 4;  // 4 electrons
    size_t K = 16;  // 16 SOs
    auto ham_par = GQCG::readFCIDUMPFile("../tests/data/lih_631g_caitlin.FCIDUMP");

    GQCG::FockSpaceProduct fock_space (16, 2);  // dim = 120
    GQCG::FCI fci (fock_space);

    // Specify solver options and solve the eigenvalue problem
    // Solve the dense FCI eigenvalue problem
    GQCG::CISolver ci_solver (fci, ham_par);
    numopt::eigenproblem::DenseSolverOptions solver_options;
    ci_solver.solve(solver_options);

    Eigen::VectorXd coef = ci_solver.get_eigenpair().get_eigenvector();
    double energy_by_eigenvalue = ci_solver.get_eigenpair().get_eigenvalue();

    // Check if the contraction energy matches the fci eigenvalue.
    GQCG::FCIRDMBuilder fci_rdm (fock_space);
    GQCG::TwoRDMs two_rdms = fci_rdm.calculate2RDMs(coef);
    GQCG::OneRDMs one_rdms = fci_rdm.calculate1RDMs(coef);

    double energy_by_contraction = ham_par.calculateEnergy(one_rdms.one_rdm, two_rdms.two_rdm);

    BOOST_CHECK(std::abs(energy_by_eigenvalue - energy_by_contraction) < 1.0e-12);
}
