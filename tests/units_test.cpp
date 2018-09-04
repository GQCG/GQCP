#define BOOST_TEST_MODULE "units"


#include "units.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise clang++ will complain



BOOST_AUTO_TEST_CASE ( bohr_angstrom ) {

    // 1 bohr = 0.529 Angstrom
    BOOST_CHECK(std::abs(GQCG::units::bohr_to_angstrom(1) - 0.529) < 1.0e-03);

    BOOST_CHECK(std::abs(GQCG::units::bohr_to_angstrom(GQCG::units::angstrom_to_bohr(0.5)) - 0.5) < 1.0e-12);
}
