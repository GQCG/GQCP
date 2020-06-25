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


#include "Utilities/type_traits.hpp"

#include <complex>


/**
 *  A header that contains special literals inside the GQCP namespace.
 */

namespace GQCP {


/**
 *  A literal for the imaginary unit.
 * 
 *  @note This feature is only added in C++14, which is why we provide it ourselves.
 */
constexpr std::complex<double> operator"" ii(long double d) {
    return std::complex<double> {0.0, static_cast<double>(d)};
}


}  // namespace GQCP
