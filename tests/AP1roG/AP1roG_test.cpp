#define BOOST_TEST_MODULE "AP1roG"


#include "AP1roG/AP1roG.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


BOOST_AUTO_TEST_CASE ( default_constructor ) {
    GQCG::AP1roG ap1rog;
}


BOOST_AUTO_TEST_CASE ( constructor ) {

    GQCG::AP1roGGeminalCoefficients g (4, 6);
    GQCG::AP1roG ap1rog (g, 0.0);
}


BOOST_AUTO_TEST_CASE ( get_geminal_coefficients ) {

    GQCG::AP1roGGeminalCoefficients g (4, 6);
    GQCG::AP1roG ap1rog (g, 0.0);
    ap1rog.get_geminal_coefficients();
}
