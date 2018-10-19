#define BOOST_TEST_MODULE "RDMCalculator_test"


#include "RDM/RDMCalculator.hpp"

#include "FockSpace/FockSpace.hpp"
#include "FockSpace/FockSpaceProduct.hpp"
#include "CISolver/CISolver.hpp"
#include "HamiltonianBuilder/DOCI.hpp"
#include "HamiltonianParameters/HamiltonianParameters_constructors.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


BOOST_AUTO_TEST_CASE ( constructor ) {

    // Test polymorphic entry for RDM (from DOCIRDMBuilder test-case).

    // Get the 1-RDM from DOCI
    size_t N = 4;  // 4 electrons
    auto ham_par = GQCG::readFCIDUMPFile("../tests/data/lih_631g_caitlin.FCIDUMP");
    size_t K = ham_par.get_K();  // 16 SO

    // Abstract pointer to test RDM
    std::shared_ptr<GQCG::BaseFockSpace> fock_space_dy(new GQCG::FockSpace(K, N/2));  // dim = 120
    GQCG::FockSpace fock_space (K, N/2);  // dim = 120

    GQCG::DOCI doci (fock_space);

    // Specify solver options and solve the eigenvalue problem
    // Solve the dense DOCI eigenvalue problem
    GQCG::CISolver ci_solver (doci, ham_par);
    numopt::eigenproblem::DenseSolverOptions solver_options;
    ci_solver.solve(solver_options);

    Eigen::VectorXd coef = ci_solver.get_eigenpair().get_eigenvector();
    
    // Check if the DOCI 1-RDM has the proper trace.
    GQCG::RDMCalculator doci_rdm (*fock_space_dy);
    GQCG::OneRDMs one_rdms = doci_rdm.calculate1RDMs(coef);

    BOOST_CHECK(std::abs(one_rdms.one_rdm.trace() - N) < 1.0e-12);
}

