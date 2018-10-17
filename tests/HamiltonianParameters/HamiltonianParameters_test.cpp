#define BOOST_TEST_MODULE "HamiltonianParameters"


#include "HamiltonianParameters/HamiltonianParameters.hpp"

#include <cpputil.hpp>

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


/*
 *  HELPER FUNCTIONS
 */
/**
 *  @return a toy 2-RDM where
 *      d(i,j,k,l) = l + 2k + 4j + 8i
 */
Eigen::Tensor<double, 4> calculateToy2RDMTensor () {
    Eigen::Tensor<double, 4> d (2, 2, 2, 2);

    for (size_t i = 0; i < 2; i++) {
        for (size_t j = 0; j < 2; j++) {
            for (size_t k = 0; k < 2; k++) {
                for (size_t l = 0; l < 2; l++) {
                    auto i_ = static_cast<double>(i);
                    auto j_ = static_cast<double>(j);
                    auto k_ = static_cast<double>(k);
                    auto l_ = static_cast<double>(l);

                    d(i,j,k,l) = l_ + 2*k_ + 4*j_ + 8*i_;
                }
            }
        }
    }

    return d;
};



/**
 *  @return toy 2-electron integrals where
 *      g(i,j,k,l) = delta_ij delta_kl - delta_il delta_jk
 */
Eigen::Tensor<double, 4> calculateToyTwoElectronIntegralsTensor () {
    Eigen::Tensor<double, 4> g (2, 2, 2, 2);
    g.setZero();

    for (size_t i = 0; i < 2; i++) {
        for (size_t j = 0; j < 2; j++) {
            for (size_t k = 0; k < 2; k++) {
                for (size_t l = 0; l < 2; l++) {
                    if ((i == j) and (k == l)) {
                        g(i,j,k,l) += 1;
                    }

                    if ((i == l) and (j == k)) {
                        g(i,j,k,l) -= 1;
                    }
                }
            }
        }
    }

    return g;
};





BOOST_AUTO_TEST_CASE ( HamiltonianParameters_constructor ) {

    // Create an AOBasis
    GQCG::Molecule water ("../tests/data/h2o.xyz");
    auto ao_basis_ptr = std::make_shared<GQCG::AOBasis>(water, "STO-3G");


    // Create One- and TwoElectronOperators (and a transformation matrix) with compatible dimensions
    size_t K = ao_basis_ptr->get_number_of_basis_functions();
    GQCG::OneElectronOperator S (Eigen::MatrixXd::Random(K, K));
    GQCG::OneElectronOperator H_core (Eigen::MatrixXd::Random(K, K));

    Eigen::Tensor<double, 4> g_tensor (K, K, K, K);
    g_tensor.setRandom();
    GQCG::TwoElectronOperator g (g_tensor);

    Eigen::MatrixXd C = Eigen::MatrixXd::Random(K, K);


    // Check if a correct constructor works
    GQCG::HamiltonianParameters random_hamiltonian_parameters (ao_basis_ptr, S, H_core, g, C);


    // Check if wrong arguments result in a throw
    GQCG::OneElectronOperator S_faulty (Eigen::MatrixXd::Random(K+1, K+1));
    GQCG::OneElectronOperator H_core_faulty (Eigen::MatrixXd::Random(K+1, K+1));

    Eigen::Tensor<double, 4> g_tensor_faulty (K+1, K+1, K+1, K+1);
    g_tensor_faulty.setRandom();
    GQCG::TwoElectronOperator g_faulty (g_tensor_faulty);

    Eigen::MatrixXd C_faulty = Eigen::MatrixXd::Random(K+1, K+1);

    BOOST_CHECK_THROW(GQCG::HamiltonianParameters (ao_basis_ptr, S_faulty, H_core, g, C), std::invalid_argument);
    BOOST_CHECK_THROW(GQCG::HamiltonianParameters (ao_basis_ptr, S, H_core_faulty, g, C), std::invalid_argument);
    BOOST_CHECK_THROW(GQCG::HamiltonianParameters (ao_basis_ptr, S, H_core, g_faulty, C), std::invalid_argument);
    BOOST_CHECK_THROW(GQCG::HamiltonianParameters (ao_basis_ptr, S, H_core, g, C_faulty), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( constructor_C ) {

    // Create dummy Hamiltonian parameters
    std::shared_ptr<GQCG::AOBasis> ao_basis;
    size_t K = 4;
    GQCG::OneElectronOperator S (Eigen::MatrixXd::Random(K, K));
    GQCG::OneElectronOperator H_core (Eigen::MatrixXd::Random(K, K));

    Eigen::Tensor<double, 4> g_tensor (K, K, K, K);
    g_tensor.setRandom();
    GQCG::TwoElectronOperator g (g_tensor);

    Eigen::MatrixXd C = Eigen::MatrixXd::Random(K, K);

    GQCG::HamiltonianParameters random_hamiltonian_parameters (ao_basis, S, H_core, g, C);


    // Check if we can create transformed Hamiltonian parameters
    Eigen::MatrixXd T = Eigen::MatrixXd::Random(K, K);
    GQCG::HamiltonianParameters transformed_random_hamiltonian_parameters (random_hamiltonian_parameters, T);
}


BOOST_AUTO_TEST_CASE ( calculate_generalized_Fock_matrix_and_super_invalid_arguments ) {

    // Initialize toy HamiltonianParameters
    std::shared_ptr<GQCG::AOBasis> ao_basis;
    GQCG::OneElectronOperator S (Eigen::MatrixXd::Zero(2, 2));
    GQCG::OneElectronOperator h (Eigen::MatrixXd::Zero(2, 2));
    Eigen::Tensor<double, 4> g_tensor (2, 2, 2, 2);
    GQCG::TwoElectronOperator g (g_tensor);
    GQCG::HamiltonianParameters ham_par (ao_basis, S, h, g, Eigen::MatrixXd::Identity(2, 2));


    // Create valid and invalid density matrices (with respect to the dimensions of the SOBasis)
    GQCG::OneRDM D_valid (Eigen::MatrixXd::Zero(2, 2));
    GQCG::OneRDM D_invalid (Eigen::MatrixXd::Zero(3, 3));

    Eigen::Tensor<double, 4> d_valid_tensor (2, 2, 2, 2);
    Eigen::Tensor<double, 4> d_invalid_tensor (3, 3, 3, 3);
    GQCG::TwoRDM d_valid (d_valid_tensor);
    GQCG::TwoRDM d_invalid (d_invalid_tensor);


    // Test a faulty function calls
    BOOST_REQUIRE_THROW(ham_par.calculateGeneralizedFockMatrix(D_invalid, d_valid), std::invalid_argument);
    BOOST_REQUIRE_THROW(ham_par.calculateGeneralizedFockMatrix(D_valid, d_invalid), std::invalid_argument);

    BOOST_REQUIRE_THROW(ham_par.calculateSuperGeneralizedFockMatrix(D_invalid, d_valid), std::invalid_argument);
    BOOST_REQUIRE_THROW(ham_par.calculateSuperGeneralizedFockMatrix(D_valid, d_invalid), std::invalid_argument);


    // Test correct function calls
    ham_par.calculateGeneralizedFockMatrix(D_valid, d_valid);
    ham_par.calculateSuperGeneralizedFockMatrix(D_valid, d_valid);
}


BOOST_AUTO_TEST_CASE ( calculate_generalized_Fock_matrix_and_super ) {

    // We test the function by a manual calculation of nonsensical toy 1- and 2-RDMS and one- and two-electron integrals
    // Set up the toy 1- and 2-RDMs
    Eigen::MatrixXd D_matrix (2, 2);
    D_matrix << 0, 1,
                2, 3;
    GQCG::OneRDM D (D_matrix);

    GQCG::TwoRDM d (calculateToy2RDMTensor());

    // Set up the toy SOBasis
    std::shared_ptr<GQCG::AOBasis> ao_basis;
    GQCG::OneElectronOperator S (Eigen::MatrixXd::Zero(2, 2));
    Eigen::MatrixXd h_matrix (2, 2);
    h_matrix << 1, 0,
                0, 1;
    GQCG::OneElectronOperator h (h_matrix);
    GQCG::TwoElectronOperator g (calculateToyTwoElectronIntegralsTensor());
    GQCG::HamiltonianParameters ham_par (ao_basis, S, h, g, Eigen::MatrixXd::Identity(2, 2));


    // Construct the reference generalized Fock matrix
    Eigen::MatrixXd F_ref = Eigen::MatrixXd::Zero(2, 2);
    for (size_t p = 0; p < 2; p++) {
        for (size_t q = 0; q < 2; q++) {
            auto p_ = static_cast<double>(p);
            auto q_ = static_cast<double>(q);

            // One-electron part is simplified by manual calculation
            F_ref(p,q) += q_ + 2*p_;

            // Two-electron part is simplified by manual calculation
            for (size_t r = 0; r < 2; r++) {
                auto r_ = static_cast<double>(r);

                F_ref(p,q) += r_ + 4*q_;
                F_ref(p,q) -= q_ + 4*r_;
            }
        }
    }


    // Construct the reference super generalized Fock matrix
    Eigen::Tensor<double, 4> W_ref (2, 2, 2, 2);
    W_ref.setZero();
    for (size_t p = 0; p < 2; p++) {
        for (size_t q = 0; q < 2; q++) {
            for (size_t r = 0; r < 2; r++) {
                for (size_t s = 0; s < 2; s++) {
                    auto q_ = static_cast<double>(q);
                    auto r_ = static_cast<double>(r);

                    if (r == q) {
                        W_ref(p,q,r,s) += F_ref(p,s);
                    }

                    // One-electron part is simplified by manual calculation
                    if (s == p) {
                        W_ref(p,q,r,s) -= q_ + 2*r_;
                    }

                    // Two-electron part is simplified by manual calculation
                    if (s == p) {
                        for (size_t t = 0; t < 2; t++) {
                            auto t_ = static_cast<double>(t);

                            W_ref(p,q,r,s) += 3*t_ - 3*q_;
                        }
                    }
                }
            }
        }
    }

    BOOST_CHECK(F_ref.isApprox(ham_par.calculateGeneralizedFockMatrix(D, d).get_matrix_representation(), 1.0e-12));
    BOOST_CHECK(cpputil::linalg::areEqual(W_ref, ham_par.calculateSuperGeneralizedFockMatrix(D, d).get_matrix_representation(), 1.0e-12));
}
