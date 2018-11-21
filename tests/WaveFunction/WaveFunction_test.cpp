// This file is part of GQCG-gqcp.
//
// Copyright (C) 2017-2018  the GQCG developers
//
// GQCG-gqcp is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GQCG-gqcp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-gqcp.  If not, see <http://www.gnu.org/licenses/>.
//
#define BOOST_TEST_MODULE "WaveFunction"

#include "WaveFunction/WaveFunction.hpp"
#include "FockSpace/FockSpace.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain



BOOST_AUTO_TEST_CASE ( shannon_entropy ) {

    // Set up a test Fock space
    GQCP::FockSpace fock_space (8, 3);  // K = 8, N = 3


    // Check the entropy of a Hartree-Fock expansion
    GQCP::WaveFunction hartree_fock_expansion (fock_space, fock_space.HartreeFockExpansion());
    BOOST_CHECK(hartree_fock_expansion.calculateShannonEntropy() < 1.0e-12);  // should be 0


    // Check the maximal entropy (corresponding to a wave function with all equal coefficients different from zero)
    GQCP::WaveFunction constant_expansion (fock_space, fock_space.constantExpansion());
    double reference_entropy = std::log2(fock_space.get_dimension());  // manual derivation
    BOOST_CHECK(std::abs(constant_expansion.calculateShannonEntropy() - reference_entropy) < 1.0e-12);
}
