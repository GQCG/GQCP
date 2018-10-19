#define BOOST_TEST_MODULE "RHF"

#include "RHF/RHF.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>


BOOST_AUTO_TEST_CASE ( HOMO_LUMO_index ) {

    // For K=7 and N=10, the index of the HOMO should be 4
    size_t K = 7;
    size_t N = 10;

    BOOST_CHECK_EQUAL(GQCP::RHFHOMOIndex(N), 4);
    BOOST_CHECK_EQUAL(GQCP::RHFLUMOIndex(K, N), 5);

    BOOST_CHECK_THROW(GQCP::RHFHOMOIndex(N+1), std::invalid_argument);
    BOOST_CHECK_THROW(GQCP::RHFLUMOIndex(K, N+1), std::invalid_argument);
}
