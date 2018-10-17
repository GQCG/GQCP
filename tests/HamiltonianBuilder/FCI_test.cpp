#define BOOST_TEST_MODULE "FCI"


#include "HamiltonianBuilder/FCI.hpp"


#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain



BOOST_AUTO_TEST_CASE ( FCI_constructor ) {
    // Create a compatible Fock space
    GQCG::FockSpaceProduct fock_space (15, 3, 3);

    // Check if a correct constructor works
    BOOST_CHECK_NO_THROW(GQCG::FCI fci (fock_space));
}


BOOST_AUTO_TEST_CASE ( FCI_public_methods ) {
    // Create an AOBasis
    GQCG::Molecule water ("../tests/data/h2o.xyz");
    auto ao_basis = std::make_shared<GQCG::AOBasis>(water, "STO-3G");


    // Create random HamiltonianParameters from One- and TwoElectronOperators (and a transformation matrix) with compatible dimensions
    size_t K = ao_basis->get_number_of_basis_functions();
    GQCG::OneElectronOperator S (Eigen::MatrixXd::Random(K, K));
    GQCG::OneElectronOperator H_core (Eigen::MatrixXd::Random(K, K));
    Eigen::Tensor<double, 4> g_tensor (K, K, K, K);
    g_tensor.setRandom();
    GQCG::TwoElectronOperator g (g_tensor);
    Eigen::MatrixXd C = Eigen::MatrixXd::Random(K, K);
    GQCG::HamiltonianParameters random_hamiltonian_parameters (ao_basis, S, H_core, g, C);

    // Create a compatible Fock space
    GQCG::FockSpaceProduct fock_space (K, 3, 3);

    // Create FCI module
    GQCG::FCI random_fci (fock_space);

    // Test the public FCI methods
    Eigen::VectorXd x = random_fci.calculateDiagonal(random_hamiltonian_parameters);
    BOOST_CHECK_NO_THROW(random_fci.constructHamiltonian(random_hamiltonian_parameters));
    BOOST_CHECK_NO_THROW(random_fci.matrixVectorProduct(random_hamiltonian_parameters, x, x));

    // Create an incompatible Fock space
    GQCG::FockSpaceProduct fock_space_i (K+1, 3, 3);

    // Create FCI module
    GQCG::FCI random_fci_i (fock_space_i);
    BOOST_CHECK_THROW(random_fci_i.constructHamiltonian(random_hamiltonian_parameters), std::invalid_argument);
    BOOST_CHECK_THROW(random_fci_i.matrixVectorProduct(random_hamiltonian_parameters, x, x), std::invalid_argument);
}
