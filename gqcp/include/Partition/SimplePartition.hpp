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


#include "Partition/PartitionTraits.hpp"
#include "Utilities/CRTP.hpp"

#include <string>
#include <vector>


namespace GQCP {


/**
 *  A type specifically designed to act as a parent class for e.g. domain-like, electron-like or orbital-like partitions in order to share common functionality.
 *
 *  @tparam _DerivedPartition              The type of the partition that derives from this class, enabling CRTP and compile-time polymorphism.
 */
template <typename _DerivedPartition>
class SimplePartition:
    public CRTP<_DerivedPartition> {
public:
    // The type of partition that derives from this class, enabling CRTP and compile-time polymorphism.
    using DerivedPartition = _DerivedPartition;
    // The type of elements that are present in the partition: domains, electrons or orbitals.
    using ElementType = typename PartitionTraits<DerivedPartition>::ElementType;

protected:
    // The partition with elements such as domains, electrons or orbitals.
    std::vector<ElementType> partition;


public:
    /*
     *  MARK: Constructors
     */

    SimplePartition() {}

    // COPY CONSTRUCTOR
    SimplePartition(const SimplePartition& partition) :
        partition {partition.partitionElements()} {}

    /**
     *  Create a partition from a vector of `element`s.
     *
     *  @param partition_elements               The `element`s that the partition contains.
     */
    SimplePartition(const std::vector<ElementType>& partition_elements) :
        partition {partition_elements} {}

    /*
     *  MARK: General info
     */

    /**
     *  @return The partition string representation.
     */
    virtual std::string asString() const = 0;

    /**
     *  @return     The dimension of this partition.
     */
    size_t dimension() const { return this->partition.size(); }

    /**
     *  @param other        The other partition.
     *
     *  @return whether this DerivedPartition is equal to the other DerivedPartition.
     */
    bool operator==(const DerivedPartition& other) const {
        bool equal = true;
        for (size_t i = 0; i < this->dimension(); i++) {
            if (this->operator()(i) != other(i)) {
                equal = false;
            }
        }
        return equal;
    }

    /**
     *  @param other        The other partition.
     *
     *  @return whether this DerivedPartition is not equal to the other DerivedPartition.
     */
    bool operator!=(const DerivedPartition& other) const {
        return !this->operator==(other);
    }

    /**
     *  @param i            The partition index.
     *
     *  @return     The i-th element of the partition.
     */
    const ElementType& operator()(const size_t i) const { return this->partition[i]; }

    /**
     *  @param os       The output stream which the partition should be concatenated to.
     *  @param partition      The partition that should be concatenated to the output stream.
     *
     *  @return The updated output stream.
     */
    friend std::ostream& operator<<(std::ostream& os, const DerivedPartition& partition) {
        return os << partition.asString();
    }

    /**
     *  @return     The elements of this partition.
     */
    const std::vector<ElementType>& partitionElements() const { return this->partition; }
};


}  // namespace GQCP
