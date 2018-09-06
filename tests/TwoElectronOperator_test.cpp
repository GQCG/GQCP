#define BOOST_TEST_MODULE "TwoElectronOperator"


#include "TwoElectronOperator.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


BOOST_AUTO_TEST_CASE ( TwoElectronOperator_constructor ) {

    // Check a correct constructor
    Eigen::Tensor<double, 4> tensor (3, 3, 3, 3);
    GQCG::TwoElectronOperator O (tensor);


    // Check a faulty constructor
    Eigen::Tensor<double, 4> tensor2 (3, 3, 3, 2);
    BOOST_CHECK_THROW(GQCG::TwoElectronOperator O2 (tensor2), std::invalid_argument); }


BOOST_AUTO_TEST_CASE ( TwoElectronOperator_getters ) {

    Eigen::Tensor<double, 4> tensor (3, 3, 3, 3);
    GQCG::TwoElectronOperator O (tensor);

    O.get_matrix_representation();
}


BOOST_AUTO_TEST_CASE ( TwoElectronOperator_transform_trivial ) {

    // Let's test a trivial transformation: i.e. with T being a unit matrix
    Eigen::Tensor<double, 4> g (3, 3, 3, 3);
    g.setRandom();
    GQCG::TwoElectronOperator G (g);

    Eigen::MatrixXd T = Eigen::MatrixXd::Identity(3, 3);
    G.transform(T);

    BOOST_CHECK(cpputil::linalg::areEqual(g, g_transformed, 1.0e-12));
}
