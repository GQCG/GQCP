#define BOOST_TEST_MODULE "CISolver"


#include "CISolver/CISolver.hpp"
#include "HamiltonianBuilder/DOCI.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


BOOST_AUTO_TEST_CASE ( Solver_constructor ) {
    // Create an AOBasis
    GQCP::Molecule water ("../tests/data/h2o.xyz");
    auto ao_basis = std::make_shared<GQCP::AOBasis>(water, "STO-3G");


    // Create random HamiltonianParameters from One- and TwoElectronOperators (and a transformation matrix) with compatible dimensions
    size_t K = ao_basis->get_number_of_basis_functions();
    GQCP::OneElectronOperator S (Eigen::MatrixXd::Random(K, K));
    GQCP::OneElectronOperator H_core (Eigen::MatrixXd::Random(K, K));
    Eigen::Tensor<double, 4> g_tensor (K, K, K, K);
    g_tensor.setRandom();
    GQCP::TwoElectronOperator g (g_tensor);
    Eigen::MatrixXd C = Eigen::MatrixXd::Random(K, K);
    GQCP::HamiltonianParameters random_hamiltonian_parameters (ao_basis, S, H_core, g, C);

    // Create a compatible Fock space
    GQCP::FockSpace fock_space (K, 3);

    // Create DOCI module
    GQCP::DOCI random_doci (fock_space);

    // Test the constructor
    BOOST_CHECK_NO_THROW(GQCP::CISolver ci_solver (random_doci, random_hamiltonian_parameters));

    // Create an incompatible Fock space
    GQCP::FockSpace fock_space_i (K+1, 3);

    // Create DOCI module
    GQCP::DOCI random_doci_i (fock_space_i);

    // Test faulty constructor
    BOOST_CHECK_THROW(GQCP::CISolver ci_solver (random_doci_i, random_hamiltonian_parameters), std::invalid_argument);
}
