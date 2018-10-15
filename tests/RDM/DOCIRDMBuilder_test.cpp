#define BOOST_TEST_MODULE "DOCI_RDM_test"



#include "RDM/DOCIRDMBuilder.hpp"

#include "CISolver/CISolver.hpp"
#include "HamiltonianBuilder/DOCI.hpp"
#include "HamiltonianParameters/HamiltonianParameters_constructors.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>


BOOST_AUTO_TEST_CASE ( lih_1RDM_trace ) {

    // Test if the trace of the 1-RDM gives N

    // Get the 1-RDM from DOCI
    size_t N = 4;  // 4 electrons
    size_t K = 16;  // 16 SOs
    auto ham_par = GQCG::readFCIDUMPFile("../tests/data/lih_631g_caitlin.FCIDUMP");

    GQCG::FockSpace fock_space (16, 2);  // dim = 120
    GQCG::DOCI doci (fock_space);

    // Specify solver options and solve the eigenvalue problem
    // Solve the dense DOCI eigenvalue problem
    GQCG::CISolver ci_solver (doci, ham_par);
    numopt::eigenproblem::DenseSolverOptions solver_options;
    ci_solver.solve(solver_options);

    Eigen::VectorXd coef = ci_solver.get_eigenpair().get_eigenvector();

    GQCG::DOCIRDMBuilder doci_rdm (fock_space);
    GQCG::OneRDM one_rdm = doci_rdm.construct1RDM(coef);

    BOOST_CHECK(std::abs(one_rdm.trace() - N) < 1.0e-12);
}


BOOST_AUTO_TEST_CASE ( lih_2RDM_trace ) {

    // Test if the trace of the 2-RDM (d_ppqq) gives N(N-1)


    // Get the 2-RDM from DOCI
    size_t N = 4;  // 4 electrons
    size_t K = 16;  // 16 SOs
    auto ham_par = GQCG::readFCIDUMPFile("../tests/data/lih_631g_caitlin.FCIDUMP");

    GQCG::FockSpace fock_space (16, 2);  // dim = 120
    GQCG::DOCI doci (fock_space);

    // Specify solver options and solve the eigenvalue problem
    // Solve the dense DOCI eigenvalue problem
    GQCG::CISolver ci_solver (doci, ham_par);
    numopt::eigenproblem::DenseSolverOptions solver_options;
    ci_solver.solve(solver_options);

    Eigen::VectorXd coef = ci_solver.get_eigenpair().get_eigenvector();

    GQCG::DOCIRDMBuilder doci_rdm (fock_space);
    GQCG::TwoRDM two_rdm = doci_rdm.construct2RDM(coef);

    BOOST_CHECK(std::abs(two_rdm.trace() - N*(N-1)) < 1.0e-12);
}


BOOST_AUTO_TEST_CASE ( lih_1RDM_2RDM_trace_DOCI ) {

    // Test if the relevant 2-RDM trace gives the 1-RDM for DOCI


    // Get the 1- and 2-RDMs from DOCI
    size_t N = 4;  // 4 electrons
    size_t K = 16;  // 16 SOs
    auto ham_par = GQCG::readFCIDUMPFile("../tests/data/lih_631g_caitlin.FCIDUMP");

    GQCG::FockSpace fock_space (16, 2);  // dim = 120
    GQCG::DOCI doci (fock_space);

    // Specify solver options and solve the eigenvalue problem
    // Solve the dense DOCI eigenvalue problem
    GQCG::CISolver ci_solver (doci, ham_par);
    numopt::eigenproblem::DenseSolverOptions solver_options;
    ci_solver.solve(solver_options);

    Eigen::VectorXd coef = ci_solver.get_eigenpair().get_eigenvector();

    GQCG::DOCIRDMBuilder doci_rdm (fock_space);
    GQCG::TwoRDM two_rdm = doci_rdm.construct2RDM(coef);
    GQCG::OneRDM one_rdm = doci_rdm.construct1RDM(coef);


    Eigen::MatrixXd D_from_reduction = (1.0/(N-1)) * two_rdm.reduce_2RDM();
    BOOST_CHECK(one_rdm.get_matrix_representation().isApprox(D_from_reduction, 1.0e-12));
}


BOOST_AUTO_TEST_CASE ( lih_energy_RDM_contraction_DOCI ) {

    // Test if the contraction of the 1- and 2-RDMs with the one- and two-electron integrals gives the DOCI energy

    size_t N = 4;  // 4 electrons
    size_t K = 16;  // 16 SOs
    auto ham_par = GQCG::readFCIDUMPFile("../tests/data/lih_631g_caitlin.FCIDUMP");

    GQCG::FockSpace fock_space (16, 2);  // dim = 120
    GQCG::DOCI doci (fock_space);

    // Specify solver options and solve the eigenvalue problem
    // Solve the dense DOCI eigenvalue problem
    GQCG::CISolver ci_solver (doci, ham_par);
    numopt::eigenproblem::DenseSolverOptions solver_options;
    ci_solver.solve(solver_options);

    Eigen::VectorXd coef = ci_solver.get_eigenpair().get_eigenvector();
    double energy_by_eigenvalue = ci_solver.get_eigenpair().get_eigenvalue();

    GQCG::DOCIRDMBuilder doci_rdm (fock_space);
    GQCG::TwoRDM two_rdm = doci_rdm.construct2RDM(coef);
    GQCG::OneRDM one_rdm = doci_rdm.construct1RDM(coef);

    Eigen::MatrixXd D = one_rdm.get_matrix_representation();
    Eigen::MatrixXd h = ham_par.get_h().get_matrix_representation();
    double energy_by_contraction = (h * D).trace();

    Eigen::Tensor<double, 4> d = two_rdm.get_tensor_representation();
    Eigen::Tensor<double, 4> g = ham_par.get_g().get_matrix_representation();

    // Specify the contractions for the relevant contraction of the two-electron integrals and the 2-RDM
    //      0.5 g(p q r s) d(p q r s)
    Eigen::array<Eigen::IndexPair<int>, 4> contractions = {Eigen::IndexPair<int>(0,0), Eigen::IndexPair<int>(1,1), Eigen::IndexPair<int>(2,2), Eigen::IndexPair<int>(3,3)};
    //      Perform the contraction
    Eigen::Tensor<double, 0> contraction = 0.5 * g.contract(d, contractions);

    // As the contraction is a scalar (a tensor of rank 0), we should access by (0).
    energy_by_contraction += contraction(0);


    BOOST_CHECK(std::abs(energy_by_eigenvalue - energy_by_contraction) < 1.0e-12);
}
