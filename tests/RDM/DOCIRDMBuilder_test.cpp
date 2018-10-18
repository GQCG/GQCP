#define BOOST_TEST_MODULE "DOCI_RDM_test"



#include "RDM/RDMBuilder.hpp"

#include "CISolver/CISolver.hpp"
#include "HamiltonianBuilder/DOCI.hpp"
#include "HamiltonianParameters/HamiltonianParameters_constructors.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>


BOOST_AUTO_TEST_CASE ( lih_1RDM_trace ) {

    // Test if the trace of the 1-RDM gives N

    // Get the 1-RDM from DOCI
    size_t N = 4;  // 4 electrons
    auto ham_par = GQCG::readFCIDUMPFile("../tests/data/lih_631g_caitlin.FCIDUMP");
    size_t K = ham_par.get_K();  // 16 SO

    GQCG::FockSpace fock_space (K, N/2);  // dim = 120
    GQCG::DOCI doci (fock_space);

    // Specify solver options and solve the eigenvalue problem
    // Solve the dense DOCI eigenvalue problem
    GQCG::CISolver ci_solver (doci, ham_par);
    numopt::eigenproblem::DenseSolverOptions solver_options;
    ci_solver.solve(solver_options);

    Eigen::VectorXd coef = ci_solver.get_eigenpair().get_eigenvector();

    // Check if the DOCI 1-RDM has the proper trace.
    GQCG::RDMBuilder doci_rdm (fock_space);
    GQCG::OneRDMs one_rdms = doci_rdm.calculate1RDMs(coef);

    BOOST_CHECK(std::abs(one_rdms.one_rdm.trace() - N) < 1.0e-12);
}


BOOST_AUTO_TEST_CASE ( lih_2RDM_trace ) {

    // Test if the trace of the 2-RDM (d_ppqq) gives N(N-1)


    // Get the 1-RDM from DOCI
    size_t N = 4;  // 4 electrons
    auto ham_par = GQCG::readFCIDUMPFile("../tests/data/lih_631g_caitlin.FCIDUMP");
    size_t K = ham_par.get_K();  // 16 SO

    GQCG::FockSpace fock_space (K, N/2);  // dim = 120
    GQCG::DOCI doci (fock_space);

    // Specify solver options and solve the eigenvalue problem
    // Solve the dense DOCI eigenvalue problem
    GQCG::CISolver ci_solver (doci, ham_par);
    numopt::eigenproblem::DenseSolverOptions solver_options;
    ci_solver.solve(solver_options);

    Eigen::VectorXd coef = ci_solver.get_eigenpair().get_eigenvector();

    // Check if the 2-RDM has the proper trace.
    GQCG::RDMBuilder doci_rdm (fock_space);
    GQCG::TwoRDMs two_rdms = doci_rdm.calculate2RDMs(coef);

    BOOST_CHECK(std::abs(two_rdms.two_rdm.trace() - N*(N-1)) < 1.0e-12);
}


BOOST_AUTO_TEST_CASE ( lih_1RDM_2RDM_trace_DOCI ) {

    // Test if the relevant 2-RDM trace gives the 1-RDM for DOCI


    // Get the 1-RDM from DOCI
    size_t N = 4;  // 4 electrons
    auto ham_par = GQCG::readFCIDUMPFile("../tests/data/lih_631g_caitlin.FCIDUMP");
    size_t K = ham_par.get_K();  // 16 SO

    GQCG::FockSpace fock_space (K, N/2);  // dim = 120
    GQCG::DOCI doci (fock_space);

    // Specify solver options and solve the eigenvalue problem
    // Solve the dense DOCI eigenvalue problem
    GQCG::CISolver ci_solver (doci, ham_par);
    numopt::eigenproblem::DenseSolverOptions solver_options;
    ci_solver.solve(solver_options);

    Eigen::VectorXd coef = ci_solver.get_eigenpair().get_eigenvector();

    // Check if the 2-RDM contraction matches the reduction.
    GQCG::RDMBuilder doci_rdm (fock_space);
    GQCG::TwoRDMs two_rdms = doci_rdm.calculate2RDMs(coef);
    GQCG::OneRDMs one_rdms = doci_rdm.calculate1RDMs(coef);


    Eigen::MatrixXd D_from_reduction = (1.0/(N-1)) * two_rdms.two_rdm.reduce();
    BOOST_CHECK(one_rdms.one_rdm.get_matrix_representation().isApprox(D_from_reduction, 1.0e-12));
}


BOOST_AUTO_TEST_CASE ( lih_energy_RDM_contraction_DOCI ) {

    // Test if the contraction of the 1- and 2-RDMs with the one- and two-electron integrals gives the DOCI energy

    // Get the 1-RDM from DOCI
    size_t N = 4;  // 4 electrons
    auto ham_par = GQCG::readFCIDUMPFile("../tests/data/lih_631g_caitlin.FCIDUMP");
    size_t K = ham_par.get_K();  // 16 SO

    GQCG::FockSpace fock_space (K, N/2);  // dim = 120
    GQCG::DOCI doci (fock_space);

    // Specify solver options and solve the eigenvalue problem
    // Solve the dense DOCI eigenvalue problem
    GQCG::CISolver ci_solver (doci, ham_par);
    numopt::eigenproblem::DenseSolverOptions solver_options;
    ci_solver.solve(solver_options);

    Eigen::VectorXd coef = ci_solver.get_eigenpair().get_eigenvector();
    double energy_by_eigenvalue = ci_solver.get_eigenpair().get_eigenvalue();

    // Check if the contraction energy matches the doci eigenvalue.
    GQCG::RDMBuilder doci_rdm (fock_space);
    GQCG::TwoRDMs two_rdms = doci_rdm.calculate2RDMs(coef);
    GQCG::OneRDMs one_rdms = doci_rdm.calculate1RDMs(coef);

    double energy_by_contraction = ham_par.calculateEnergy(one_rdms.one_rdm, two_rdms.two_rdm);

    BOOST_CHECK(std::abs(energy_by_eigenvalue - energy_by_contraction) < 1.0e-12);
}
