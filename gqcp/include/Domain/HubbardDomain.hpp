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


#include "Domain/DiscreteDomain.hpp"
#include "QuantumChemical/SpinResolved.hpp"


namespace GQCP {


/**
 * A Hubbard domain as a collection of sites.
 * The sites {`i'} that are present in the domain are represented by a set bit at the corresponding indices `i'.
 */
class HubbardDomain:
    public DiscreteDomain {
public:
    using DiscreteDomain::DiscreteDomain;

    /**
     * Construct the restricted projection matrix of this Hubbard Domain.
     *
     * @param L     Number of sites in the Hubbard system.
     *
     * @return The projection matrix for this Hubbard domain.
     */
    RTransformation<double> RProjectionMatrix(const size_t& L) const {
        return RTransformation<double> {SquareMatrix<double>::PartitionMatrix(this->domainIndices(), L)};
    }

    /**
     * Construct the unrestricted projection matrix of this Hubbard Domain.
     *
     * @param L     Number of sites in the Hubbard system.
     *
     * @return The projection matrix for this Hubbard domain.
     */
    UTransformation<double> UProjectionMatrix(const size_t& L) const {
        return UTransformation<double> {UTransformationComponent<double> {SquareMatrix<double>::PartitionMatrix(this->domainIndices(), L)}, UTransformationComponent<double> {SquareMatrix<double>::PartitionMatrix(this->domainIndices(), L)}};
    }
};


}  // namespace GQCP
