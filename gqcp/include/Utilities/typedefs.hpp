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


#include <complex>
#include <cstdlib>
#include <type_traits>
#include <vector>


/**
 *  A header that contains general typedefs inside the GQCP namespace
 */


namespace GQCP {


/*
 *  SCALARS
 */
using cd = std::complex<double>;


/*
 *  TEMPLATE ALIASES
 */
template <typename T, typename U>
using sum_t = decltype(std::declval<T>() + std::declval<U>());

template <typename T, typename U>
using product_t = decltype(std::declval<T>() * std::declval<U>());


}  // namespace GQCP
