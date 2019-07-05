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
// 
#ifndef GQCP_LIBINTONEELECTRONINTEGRALBUFFER_HPP
#define GQCP_LIBINTONEELECTRONINTEGRALBUFFER_HPP


#include "Basis/BaseOneElectronIntegralBuffer.hpp"


namespace GQCP {


/**
 *  A buffer for storing libint one-electron integrals
 * 
 *  @tparam _N              the number of components the operator has
 */
template <size_t _N>
class LibintOneElectronIntegralBuffer : public BaseOneElectronIntegralBuffer<double, _N> {
public:
    using Scalar = double;  // the scalar representation of an integral for libint is always a real number
    static constexpr auto N = _N;  // the number of components the operator has
};


}  // namespace GQCP



#endif  // GQCP_LIBINTONEELECTRONINTEGRALBUFFER_HPP
