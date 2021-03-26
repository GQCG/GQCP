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


#include "Mathematical/Functions/EvaluableLinearCombination.hpp"
#include "QuantumChemical/SpinResolvedBase.hpp"


namespace GQCP {


/**
 *  @tparam _Scalar                 The scalar representation of one of the expansion coefficients.
 *  @tparam _BasisFunction          The type of the underlying basis functions.
 */
template <typename _Scalar, typename _BasisFunction>
class Spinor:
    public SpinResolvedBase<EvaluableLinearCombination<_Scalar, _BasisFunction>, Spinor<_Scalar, _BasisFunction>> {
public:
    // The scalar representation of one of the expansion coefficients.
    using Scalar = _Scalar;

    // The type of the underlying basis functions.
    using BasisFunction = _BasisFunction;
};


}  // namespace GQCP
