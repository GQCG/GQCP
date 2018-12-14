// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2018  the GQCG developers
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
#ifndef GQCP_UNITS_HPP
#define GQCP_UNITS_HPP


namespace GQCP {
namespace units {


namespace constants {

struct CODATA2014 {
    constexpr static double angstrom_per_bohr = 0.52917721067;
    constexpr static double bohr_per_angstrom = 1 / angstrom_per_bohr;
};

}  // namespace constants



/**
 *  @param value_in_bohr        a distance expressed in bohr
 *
 *  @return the value in angstrom
 */
inline double bohr_to_angstrom(double value_in_bohr) { return value_in_bohr * constants::CODATA2014::angstrom_per_bohr; }

/**
 *  @param value_in_angstrom        a distance expressed in angstrom
 *
 *  @return the value in bohr (a.u.)
 */
inline double angstrom_to_bohr(double value_in_angstrom) { return value_in_angstrom * constants::CODATA2014::bohr_per_angstrom; }



}  // namespace units
}  // namespace GQCP


#endif  // GQCP_UNITS_HPP
