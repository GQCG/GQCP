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


/**
 *  A class that represents a spinor basis without any restrictions on the expansion of the alpha and beta components in terms of the underlying (possibly different) scalar bases
 * 
 *  @tparam _Shell                      the type of shell the underlying scalar bases contain
 *  @tparam _TransformationScalar       the scalar type of the expansion coefficients
 */
template <typename _TransformationScalar, typename _Shell>
class RSpinorBasis {