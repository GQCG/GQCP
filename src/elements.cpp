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
#include "elements.hpp"

#include <boost/bimap.hpp>


namespace GQCP {
namespace elements {


/*
 *  STATIC VARIABLES
 */

/**
 *  Create a vector of bimap values, which provides .begin() and .end() iterators
 *
 *  Adapted from https://stackoverflow.com/a/20290421/7930415
 */
std::vector<boost::bimap<std::string, size_t>::value_type> elements_list {
    {"H",  1},
    {"He", 2},
    {"Li", 3},
    {"Be", 4},
    {"B",  5},
    {"C",  6},
    {"N",  7},
    {"O",  8},
    {"F",  9},
    {"Ne", 10},
    {"Na", 11},
    {"Mg", 12},
    {"Al", 13},
    {"Si", 14},
    {"P",  15},
    {"S",  16},
    {"Cl", 17},
    {"Ar", 18}
};


/**
 *  Use the .begin() and .end() iterators of a std::vector to construct a boost::bimap
 *
 *  Adapted from https://stackoverflow.com/a/20290421/7930415
 */
static const boost::bimap<std::string, size_t> periodic_table (elements_list.begin(), elements_list.end());



/*
 * FUNCTION IMPLEMENTATIONS
 */

/**
 *  Given a @param symbol for the name of an element, @return its atomic number
 */
size_t elementToAtomicNumber(const std::string& symbol) {

    return periodic_table.left.at(symbol);
}


/**
 *  Given an @param atomic_number, @return the symbol of the corresponding element
 */
const std::string& atomicNumberToElement(size_t atomic_number) {

    return periodic_table.right.at(atomic_number);
}


}  // namespace elements
}  // namespace GQCP
