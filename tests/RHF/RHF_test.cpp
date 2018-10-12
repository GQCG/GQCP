#define BOOST_TEST_MODULE "RHF"

#include "RHF/RHF.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>


BOOST_AUTO_TEST_CASE ( homo_lumo ) {

    // For K=7 and N=10, the index of the HOMO should be 4
    size_t K = 7;
    size_t N = 10;

    BOOST_CHECK_EQUAL(GQCG::HOMOIndex(N), 4);
    BOOST_CHECK_EQUAL(rhf.LUMOIndex(K, N), 5);
}
