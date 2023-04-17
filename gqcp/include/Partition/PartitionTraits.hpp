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


#include <cstdlib>


namespace GQCP {


/**
 *  A type that provides compile-time information on the type of elements in a partition that is otherwise not accessible through a public class alias.
 */
template <typename DerivedPartition>
struct PartitionTraits {};

template <typename DomainType>
class DomainPartition;

template <typename DomainType>
struct PartitionTraits<DomainPartition<DomainType>> {
    // The type of elements that are present in the domain.
    using ElementType = DomainType;
};


}  // namespace GQCP
