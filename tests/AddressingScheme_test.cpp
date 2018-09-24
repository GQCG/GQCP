#define BOOST_TEST_MODULE "AddressingScheme"


#include "AddressingScheme.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>



BOOST_AUTO_TEST_CASE ( vertex_weights_K5_N3 ) {

    // Let's test an addressing scheme for K=5 and N=3 (5 MOs and 3 alpha electrons)
    GQCG::AddressingScheme addressing_scheme = GQCG::AddressingScheme(5, 3);

    GQCG::Matrixu ref_vertex_weights = {{1, 0, 0, 0},
                                        {1, 1, 0, 0},
                                        {1, 2, 1, 0},
                                        {0, 3, 3, 1},
                                        {0, 0, 6, 4},
                                        {0, 0, 0, 10}};

    BOOST_CHECK(ref_vertex_weights == addressing_scheme.get_vertex_weights());
}
