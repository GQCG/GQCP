#define BOOST_TEST_MODULE "OneElectronOperator"


#include "OneElectronOperator.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


BOOST_AUTO_TEST_CASE ( OneElectronOperator_constructor ) {

    // Check a correct constructor
    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(4, 4);
    GQCG::OneElectronOperator O (matrix);


    // Check a faulty constructor
    Eigen::MatrixXd matrix2 = Eigen::MatrixXd::Zero(3, 4);
    BOOST_CHECK_THROW(GQCG::OneElectronOperator O2 (matrix2), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( OneElectronOperator_getters ) {

    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(4, 4);
    GQCG::OneElectronOperator O (matrix);

    O.get_matrix_representation();
}
