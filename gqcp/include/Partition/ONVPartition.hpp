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

#pragma once


#include "ONVBasis/SpinResolvedONV.hpp"
#include "ONVBasis/SpinUnresolvedONV.hpp"
#include "Partition/DiscreteDomainPartition.hpp"
#include "Partition/SimplePartition.hpp"

#include <bitset>

#include <cmath>

namespace GQCP {


/**
 * An ONV partitioned according to a set of discrete domains, i.e. a `DiscreteDomainPartition`.
 */
template <typename ONV>
class ONVPartition:
    public SimplePartition<ONVPartition<ONV>> {

private:
    int phase_factor;  // Phase factor arising from the anti-commutation rules when reordering the creation operators by domain.

public:
    /*
     *  MARK: Constructors
     */

    // Inherit `SimplePartition`'s constructors.
    using SimplePartition<ONVPartition<ONV>>::SimplePartition;

    /**
     *  Partition a spin-unresolved ONV according to a discrete domain partition.
     *
     *  @param domain_partition          The discrete domain partition.
     *  @param onv                       The spin-unresolved ONV.
     */
    ONVPartition(const DiscreteDomainPartition& domain_partition, const SpinUnresolvedONV& onv) {
        const auto K = onv.numberOfSpinors();
        const auto D = domain_partition.dimension();

        // Check whether the ONV has the same dimension as each domain in the discrete domain partition.                                                                                                                                                                          discrete domain.
        if (domain_partition(0).dimension() != K) {
            throw std::invalid_argument("ONVPartition::ONVPartition(const DiscreteDomainPartition&, const GQCP::SpinUnresolvedONV&): Please provide an ONV that has the same dimension as each domain within the domain partition.");
        }

        this->partition.reserve(D);
        int phase_factor = 1;
        // Create a domain-related ONV for each domain. The dimension of the ONV is equal to the number of elements in the domain.
        for (size_t d = 0; d < D; d++) {
            const auto domain_indices = domain_partition(d).domainIndices();

            size_t unsigned_representation = 0UL;
            size_t p_domain = 0;
            // Loop over all elements in the domain. If the domain element at index `i` contains an electron, keep track of the anticommutation rules.
            for (const size_t& p : domain_indices) {
                if (onv.isOccupied(p)) {
                    // Adjust unsigned representation of the domain ONV by setting a bit at the correct position.
                    unsigned_representation ^= 1UL << p_domain++;
                    // Check anti-commutation rules such that the creation operators are reordered according to the domains.
                    for (size_t other_d = d + 1; other_d < D; other_d++) {
                        auto occupied_in_other_d = domain_partition(other_d).unsignedRepresentation() & onv.unsignedRepresentation();
                        auto occupied_before_p = occupied_in_other_d & ~(-(1UL << p));
                        phase_factor *= std::pow(-1, __builtin_popcountll(occupied_before_p) % 2);
                    }
                } else {
                    // If the domain index at `i` is not occupied with an electron, add a zero to the unsigned representation and move on.
                    ++p_domain;
                }
            }

            this->partition.push_back(SpinUnresolvedONV(domain_indices.size(), __builtin_popcountll(unsigned_representation), unsigned_representation));
        }
        this->phase_factor = phase_factor;
    }

    /**
     *  Partition a spin-resolved ONV according to a discrete domain partition.
     *
     *  @param domain_partition          The discrete domain partition.
     *  @param onv                       The spin-resolved ONV.
     */
    ONVPartition(const DiscreteDomainPartition& domain_partition, const SpinResolvedONV& onv) {
        const auto onv_partition_alpha = ONVPartition<SpinUnresolvedONV>(domain_partition, onv.onv(Spin::alpha));
        const auto onv_partition_beta = ONVPartition<SpinUnresolvedONV>(domain_partition, onv.onv(Spin::beta));

        const auto D = domain_partition.dimension();

        this->partition.reserve(D);
        // Each domain in the domain partition will have a spin-resolved ONV.
        for (size_t d = 0; d < D; d++) {
            this->partition.push_back(GQCP::SpinResolvedONV(onv_partition_alpha(d), onv_partition_beta(d)));
        }
        this->phase_factor = onv_partition_alpha.phaseFactor() * onv_partition_beta.phaseFactor();  // Phase factor is the product of the phase factor of each spin-component separately.
    }

    /*
     *  MARK: General info
     */

    /**
     *  @return The ONV partition string representation.
     */
    std::string asString() const {
        std::string s;
        for (size_t i = 0; i < this->dimension(); i++) {
            s += this->operator()(i).asString() + " / ";
        }
        return s.substr(0, s.length() - 3);
    }

    /**
     *  @return The phase factor of the ONV partition.
     */
    int phaseFactor() const { return this->phase_factor; }
};


}  // namespace GQCP