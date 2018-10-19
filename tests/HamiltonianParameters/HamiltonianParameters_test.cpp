#define BOOST_TEST_MODULE "HamiltonianParameters"


#include "HamiltonianParameters/HamiltonianParameters.hpp"


#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain



BOOST_AUTO_TEST_CASE ( HamiltonianParameters_constructor ) {

    // Create an AOBasis
    GQCP::Molecule water ("../tests/data/h2o.xyz");
    auto ao_basis_ptr = std::make_shared<GQCP::AOBasis>(water, "STO-3G");


    // Create One- and TwoElectronOperators (and a transformation matrix) with compatible dimensions
    size_t K = ao_basis_ptr->get_number_of_basis_functions();
    GQCP::OneElectronOperator S (Eigen::MatrixXd::Random(K, K));
    GQCP::OneElectronOperator H_core (Eigen::MatrixXd::Random(K, K));

    Eigen::Tensor<double, 4> g_tensor (K, K, K, K);
    g_tensor.setRandom();
    GQCP::TwoElectronOperator g (g_tensor);

    Eigen::MatrixXd C = Eigen::MatrixXd::Random(K, K);


    // Check if a correct constructor works
    GQCP::HamiltonianParameters random_hamiltonian_parameters (ao_basis_ptr, S, H_core, g, C);


    // Check if wrong arguments result in a throw
    GQCP::OneElectronOperator S_faulty (Eigen::MatrixXd::Random(K+1, K+1));
    GQCP::OneElectronOperator H_core_faulty (Eigen::MatrixXd::Random(K+1, K+1));

    Eigen::Tensor<double, 4> g_tensor_faulty (K+1, K+1, K+1, K+1);
    g_tensor_faulty.setRandom();
    GQCP::TwoElectronOperator g_faulty (g_tensor_faulty);

    Eigen::MatrixXd C_faulty = Eigen::MatrixXd::Random(K+1, K+1);

    BOOST_CHECK_THROW(GQCP::HamiltonianParameters (ao_basis_ptr, S_faulty, H_core, g, C), std::invalid_argument);
    BOOST_CHECK_THROW(GQCP::HamiltonianParameters (ao_basis_ptr, S, H_core_faulty, g, C), std::invalid_argument);
    BOOST_CHECK_THROW(GQCP::HamiltonianParameters (ao_basis_ptr, S, H_core, g_faulty, C), std::invalid_argument);
    BOOST_CHECK_THROW(GQCP::HamiltonianParameters (ao_basis_ptr, S, H_core, g, C_faulty), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( constructor_C ) {

    // Create dummy Hamiltonian parameters
    std::shared_ptr<GQCP::AOBasis> ao_basis;
    size_t K = 4;
    GQCP::OneElectronOperator S (Eigen::MatrixXd::Random(K, K));
    GQCP::OneElectronOperator H_core (Eigen::MatrixXd::Random(K, K));

    Eigen::Tensor<double, 4> g_tensor (K, K, K, K);
    g_tensor.setRandom();
    GQCP::TwoElectronOperator g (g_tensor);

    Eigen::MatrixXd C = Eigen::MatrixXd::Random(K, K);

    GQCP::HamiltonianParameters random_hamiltonian_parameters (ao_basis, S, H_core, g, C);


    // Check if we can create transformed Hamiltonian parameters
    Eigen::MatrixXd T = Eigen::MatrixXd::Random(K, K);
    GQCP::HamiltonianParameters transformed_random_hamiltonian_parameters (random_hamiltonian_parameters, T);
}
