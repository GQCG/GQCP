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
 *  Extensions of the <type_traits> header that haven't been included in C++11
 */


#include <type_traits>



/**
 *  An implementation of enable_if_t.
 */

namespace GQCP {


template <bool B, typename T = void>
using enable_if_t = typename std::enable_if<B, T>::type;  // only in C++14


template <class T, class U>
constexpr bool is_same_v = std::is_same<T, U>::value;  // only in C++17

}  // namespace GQCP
