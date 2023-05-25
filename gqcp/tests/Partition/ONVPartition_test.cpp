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

#define BOOST_TEST_MODULE "ONVPartition"

#include <boost/test/unit_test.hpp>

#include "Partition/ONVPartition.hpp"

/**
 *  Test the spin-unresolved ONVs when partitioned according to a discrete domain partitions.
 */
BOOST_AUTO_TEST_CASE(spin_unresolved) {

    // 001011 | 010000 | 100100
    const GQCP::DiscreteDomainPartition domain_partition {std::vector<size_t> {11, 16, 36}, 6};
    // 011101
    const GQCP::SpinUnresolvedONV onv {6, 4, 29};

    const GQCP::ONVPartition<GQCP::SpinUnresolvedONV> onv_partition {domain_partition, onv};
    const std::vector<std::string> ref_onvs = {"101", "1", "10"};

    for (size_t i = 0; i < onv_partition.dimension(); ++i) {
        BOOST_CHECK_EQUAL(onv_partition(i).asString(), ref_onvs[i]);
    }
    BOOST_CHECK_EQUAL(onv_partition.phaseFactor(), 1);

    // 0100000001010001 | 1000010000000010 | 0000000100000100 | 0011100000001000 | 0000001010100000
    const GQCP::DiscreteDomainPartition other_domain_partition {std::vector<size_t> {16465, 33794, 260, 14344, 672}, 16};
    const auto other_onv = GQCP::SpinUnresolvedONV::FromString("0101111101001101");

    const GQCP::ONVPartition<GQCP::SpinUnresolvedONV> other_onv_partition {other_domain_partition, other_onv};
    const std::vector<std::string> ref_other_onvs = {"1011", "010", "11", "1110", "001"};

    for (size_t i = 0; i < other_onv_partition.dimension(); ++i) {
        BOOST_CHECK_EQUAL(other_onv_partition(i).asString(), ref_other_onvs[i]);
    }
    BOOST_CHECK_EQUAL(other_onv_partition.phaseFactor(), 1);
}


/**
 *  Test the spin-resolved ONVs when partitioned according to a discrete domain partitions.
 */
BOOST_AUTO_TEST_CASE(spin_resolved) {

    // 001011 | 010000 | 100100
    const GQCP::DiscreteDomainPartition domain_partition {std::vector<size_t> {11, 16, 36}, 6};
    // 011101 | 101100
    const GQCP::SpinResolvedONV onv {GQCP::SpinUnresolvedONV {6, 4, 29}, GQCP::SpinUnresolvedONV {6, 3, 44}};

    const GQCP::ONVPartition<GQCP::SpinResolvedONV> onv_partition {domain_partition, onv};
    const std::vector<std::string> ref_onvs_alpha = {"101", "1", "10"};
    const std::vector<std::string> ref_onvs_beta = {"001", "0", "11"};

    for (size_t i = 0; i < onv_partition.dimension(); ++i) {
        BOOST_CHECK_EQUAL(onv_partition(i).onv(GQCP::Spin::alpha).asString(), ref_onvs_alpha[i]);
        BOOST_CHECK_EQUAL(onv_partition(i).onv(GQCP::Spin::beta).asString(), ref_onvs_beta[i]);
    }
    BOOST_CHECK_EQUAL(onv_partition.phaseFactor(), -1);

    // 0100000001010001 | 1000010000000010 | 0000000100000100 | 0011100000001000 | 0000001010100000
    const GQCP::DiscreteDomainPartition other_domain_partition {std::vector<size_t> {16465, 33794, 260, 14344, 672}, 16};
    // 0101111101001101 | 0100000111011011
    const GQCP::SpinResolvedONV other_onv {GQCP::SpinUnresolvedONV {16, 10, 24397}, GQCP::SpinUnresolvedONV {16, 8, 16859}};

    const GQCP::ONVPartition<GQCP::SpinResolvedONV> other_onv_partition {other_domain_partition, other_onv};
    const std::vector<std::string> ref_other_onvs_alpha = {"1011", "010", "11", "1110", "001"};
    const std::vector<std::string> ref_other_onvs_beta = {"1111", "100", "01", "1000", "010"};

    for (size_t i = 0; i < other_onv_partition.dimension(); ++i) {
        BOOST_CHECK_EQUAL(other_onv_partition(i).onv(GQCP::Spin::alpha).asString(), ref_other_onvs_alpha[i]);
        BOOST_CHECK_EQUAL(other_onv_partition(i).onv(GQCP::Spin::beta).asString(), ref_other_onvs_beta[i]);
    }
    BOOST_CHECK_EQUAL(other_onv_partition.phaseFactor(), 1);
}
