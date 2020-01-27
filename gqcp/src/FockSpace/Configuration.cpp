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


#include "FockSpace/Configuration.hpp"


namespace GQCP {
    

/**
 *  @param delimiter        the delimiter between bits
 * 
 *  @return a compact representation of the configuration as string:
 *   e.g: 1100 | 1010 => 2110
 */
std::string Configuration::asCompactString(const std::string& delimiter) const {

    std::string beta = this->onv_beta.asString();
    std::string alpha = this->onv_alpha.asString();
    std::string result = std::to_string((int)(alpha[0]) + (int)(beta[0]) - 96);  // ASCII Code for 0 is 48, adding two ASCII CODES and substracting 98 will result in character integer as integer value

    for (size_t i = 1; i < beta.size(); i++) {
        result += delimiter;
        result += std::to_string((int)(alpha[i]) + (int)(beta[i]) - 96);;
    } 

    return result;
}



}  // namespace GQCP
