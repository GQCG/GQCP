#define BOOST_TEST_MODULE "AP1roGGeminalCoefficients"


#include "AP1roGGeminalCoefficients.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


BOOST_AUTO_TEST_CASE ( constructor ) {

    // Check a correct constructor
    GQCG::AP1roGGeminalCoefficients g (4, 6);

    // We can't create 4 geminals in 4 orbitals
    BOOST_CHECK_THROW(GQCG::AP1roGGeminalCoefficients(4, 4), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( constructor_vector ) {

    // Check a correct constructor
    Eigen::VectorXd g = Eigen::VectorXd::Zero(6);
    BOOST_CHECK_NO_THROW(GQCG::AP1roGGeminalCoefficients(g, 2, 5));

    // Check wrong parameters N_P and K
    BOOST_CHECK_THROW(GQCG::AP1roGGeminalCoefficients(g, 2, 4), std::invalid_argument);
    BOOST_CHECK_THROW(GQCG::AP1roGGeminalCoefficients(g, 1, 5), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( asVector ) {

    GQCG::AP1roGGeminalCoefficients g (4, 6);
    g.asVector();
}


BOOST_AUTO_TEST_CASE ( asMatrix ) {

    // For N_P=2 and K=5, we have an AP1roG geminal coefficient matrix that looks like the following matrix:
    Eigen::MatrixXd G (2, 5);
    G << 1, 0,  1, 2, 3,
         0, 1,  4, 5, 6;

    // The geminal coefficients, arranged in a vector are then represented by the following vector:
    Eigen::VectorXd g (6);
    g << 1, 2, 3, 4, 5, 6;


    GQCG::AP1roGGeminalCoefficients gem_coeff (g, 2, 5);
    BOOST_CHECK(gem_coeff.asMatrix().isApprox(G));
}
