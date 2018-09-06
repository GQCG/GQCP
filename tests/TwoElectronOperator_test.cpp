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
