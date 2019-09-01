// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2019  the GQCG developers
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
#pragma once


#include "Basis/SPBasis.hpp"
#include "Operator/SecondQuantized/SQOneElectronOperator.hpp"

#include <cstddef>


namespace GQCP {


/**
 *  A class that represents the overlap operator
 */
class OverlapOperator {
public:
    using Scalar = double;  // the scalar representation of the operator
    static constexpr size_t Components = 1;  // the number of components the operator has


public:
    // PUBLIC METHODS

    /**
     *  @tparam TransformationScalar        the scalar type of the transformation matrix that connects the scalar basis with the current single-particle 'orbitals'
     *  @tparam ShellType                   the type of shell that this scalar basis contains
     * 
     *  @param sp_basis                     the single-particle basis in which the integrals of the second-quantized operator should be expressed
     * 
     *  @return the second-quantized operator that corresponds to this first-quantized operator
     */
    template <typename TransformationScalar, typename ShellType>
    auto quantize(const SPBasis<TransformationScalar, ShellType>& sp_basis) const -> ScalarSQOneElectronOperator<product_t<Scalar, TransformationScalar>> {

        using ResultScalar = product_t<Scalar, TransformationScalar>;

        ScalarSQOneElectronOperator<ResultScalar> op ({sp_basis.scalarBasis().calculateLibintOverlapIntegrals()});  // op for 'operator'
        op.transform(sp_basis.transformationMatrix());
        return op;
    }
};


}  // namespace GQCP
