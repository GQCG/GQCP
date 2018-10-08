#define BOOST_TEST_MODULE "FockSpace"


#include "FockSpace/FockSpace.hpp"


#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


BOOST_AUTO_TEST_CASE ( FockSpace_constructor ) {

    BOOST_CHECK_NO_THROW(GQCG::FockSpace (10,5));
}


BOOST_AUTO_TEST_CASE ( FockSpace_dimension) {

    BOOST_CHECK_EQUAL(GQCG::FockSpace::calculateDimension(10, 1), 10);
    BOOST_CHECK_EQUAL(GQCG::FockSpace::calculateDimension(6, 2), 15);
    BOOST_CHECK_EQUAL(GQCG::FockSpace::calculateDimension(8, 3), 56);
}


BOOST_AUTO_TEST_CASE ( vertex_weights_K5_N3 ) {

    // Let's test an addressing scheme for K=5 and N=3 (5 MOs and 3 alpha electrons)
    GQCG::FockSpace fock_space = GQCG::FockSpace(5, 3);

    GQCG::Matrixu ref_vertex_weights = {{1, 0, 0, 0},
                                        {1, 1, 0, 0},
                                        {1, 2, 1, 0},
                                        {0, 3, 3, 1},
                                        {0, 0, 6, 4},
                                        {0, 0, 0, 10}};
    BOOST_CHECK(ref_vertex_weights == fock_space.get_vertex_weights());
}
