// This file is part of GQCG-GQCP.
//
// Copyright (C) 2017-2020  the GQCG developers
//
// GQCG-GQCP is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GQCG-GQCP is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-GQCP.  If not, see <http://www.gnu.org/licenses/>.

#define BOOST_TEST_MODULE "DiscreteDomainPartition"

#include <boost/test/unit_test.hpp>

#include "Partition/DiscreteDomainPartition.hpp"

/**
 * Test if the domains are correctly displayed as strings and vectors.
 */
BOOST_AUTO_TEST_CASE(constructor) {
    // 10100001 | 01010010 | 00001100
    const GQCP::DiscreteDomainPartition domain_partition {std::vector<size_t> {161, 82, 12}, 8};
    BOOST_CHECK_EQUAL(domain_partition.asString(), std::string("0-1-2-2-1-0-1-0"));

    const auto domain_partition_vector = domain_partition.asVector();
    const std::vector<size_t> ref_vector = {0, 1, 2, 2, 1, 0, 1, 0};
    BOOST_CHECK_EQUAL_COLLECTIONS(domain_partition_vector.begin(), domain_partition_vector.end(), ref_vector.begin(), ref_vector.end());

    BOOST_CHECK_THROW(GQCP::DiscreteDomainPartition fuzzy_domain_partition(std::vector<size_t> {163, 82, 12}, 8), std::invalid_argument);
}


/**
 * Test if the discrete domain partition is correctly constructed from its vector representation.
 */
BOOST_AUTO_TEST_CASE(vector_constructor) {
    // 10100001 | 01010010 | 00001100
    const GQCP::DiscreteDomainPartition domain_partition {std::vector<size_t> {0, 1, 2, 2, 1, 0, 1, 0}};
    BOOST_CHECK_EQUAL(domain_partition.asString(), std::string("0-1-2-2-1-0-1-0"));

    const auto domain_partition_vector = domain_partition.asVector();
    const std::vector<size_t> ref_vector = {0, 1, 2, 2, 1, 0, 1, 0};
    BOOST_CHECK_EQUAL_COLLECTIONS(domain_partition_vector.begin(), domain_partition_vector.end(), ref_vector.begin(), ref_vector.end());
}


/**
 * Test if the calculated overlap with a spin-unresolved ONV gives the correct electron partition.
 */
BOOST_AUTO_TEST_CASE(overlap_with_SpinUnresolvedONV) {
    // 10100001 | 01010010 | 00001100
    const GQCP::DiscreteDomainPartition domain_partition {std::vector<size_t> {161, 82, 12}, 8};
    // 10010110
    const GQCP::SpinUnresolvedONV onv {8, 4, 150};

    const auto overlap_electron_partition = domain_partition.overlapWithONV(onv);
    const std::vector<size_t> ref_electron_partition = {1, 2, 1};
    BOOST_CHECK_EQUAL_COLLECTIONS(overlap_electron_partition.partitionElements().begin(), overlap_electron_partition.partitionElements().end(), ref_electron_partition.begin(), ref_electron_partition.end());
}


/**
 * Test if the calculated overlap with a spin-resolved ONV gives the correct electron partition.
 */
BOOST_AUTO_TEST_CASE(overlap_with_SpinResolvedONV) {
    // 10100001 | 01010010 | 00001100
    const GQCP::DiscreteDomainPartition domain_partition {std::vector<size_t> {161, 82, 12}, 8};
    // 10010110 | 01001101
    const GQCP::SpinResolvedONV onv {GQCP::SpinUnresolvedONV {8, 4, 150}, GQCP::SpinUnresolvedONV {8, 4, 77}};

    const auto overlap_electron_partition = domain_partition.overlapWithONV(onv);
    const std::vector<size_t> ref_electron_partition_alpha = {1, 2, 1};
    const std::vector<size_t> ref_electron_partition_beta = {1, 1, 2};
    BOOST_CHECK_EQUAL_COLLECTIONS(overlap_electron_partition.alpha().partitionElements().begin(), overlap_electron_partition.alpha().partitionElements().end(), ref_electron_partition_alpha.begin(), ref_electron_partition_alpha.end());
    BOOST_CHECK_EQUAL_COLLECTIONS(overlap_electron_partition.beta().partitionElements().begin(), overlap_electron_partition.beta().partitionElements().end(), ref_electron_partition_beta.begin(), ref_electron_partition_beta.end());
}
