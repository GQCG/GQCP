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


BOOST_AUTO_TEST_CASE ( OneElectronOperator_transform_trivial ) {

    // Let's test a trivial transformation: i.e. with T being a unit matrix
    Eigen::MatrixXd h = Eigen::MatrixXd::Random(3, 3);
    GQCG::OneElectronOperator H (h);

    Eigen::MatrixXd T = Eigen::MatrixXd::Identity(3, 3);
    H.transform(T);

    BOOST_CHECK(H.get_matrix_representation().isApprox(h, 1.0e-12));
}


BOOST_AUTO_TEST_CASE ( OneElectronOperator_transform_and_inverse ) {

    // Let's test if, if we transform h with T and then with T_inverse, we get effectively do nothing
    Eigen::MatrixXd h = Eigen::MatrixXd::Random(3, 3);
    GQCG::OneElectronOperator H (h);

    Eigen::MatrixXd T (3, 3);
    T << 1,  0,  0,
         0, -2,  0,
         0,  0,  3;
    Eigen::MatrixXd T_inverse = T.inverse();


    H.transform(T);
    H.transform(T_inverse);

    BOOST_CHECK(H.get_matrix_representation().isApprox(h, 1.0e-12));
}
