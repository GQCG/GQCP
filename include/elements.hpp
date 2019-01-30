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
#ifndef GQCP_ELEMENTS_HPP
#define GQCP_ELEMENTS_HPP


#include <string>


namespace GQCP {
namespace elements {


/**
 *  @param symbol       the name of an element
 *
 *  @return the atomic number of the corresponding element
 */
size_t elementToAtomicNumber(const std::string& symbol);


/**
 *  @param atomic_number    the atomic number of an element
 *
 *  @return the symbol of the corresponding element
 */
const std::string& atomicNumberToElement(size_t atomic_number);



}  // namespace elements
}  // namespace GQCP


#endif  // GQCP_ELEMENTS_HPP
