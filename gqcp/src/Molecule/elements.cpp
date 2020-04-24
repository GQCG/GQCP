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

#include "Molecule/elements.hpp"

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
// clang-format off
std::vector<boost::bimap<std::string, size_t>::value_type> elements_list {
    {"H",     1},
    {"He",    2},
    {"Li",    3},
    {"Be",    4},
    {"B",     5},
    {"C",     6},
    {"N",     7},
    {"O",     8},
    {"F",     9},
    {"Ne",   10},
    {"Na",   11},
    {"Mg",   12},
    {"Al",   13},
    {"Si",   14},
    {"P",    15},
    {"S",    16},
    {"Cl",   17},
    {"Ar",   18},
    {"K",    19},
    {"Ca",   20},
    {"Sc",   21},
    {"Ti",   22},
    {"V",    23},
    {"Cr",   24},
    {"Mn",   25},
    {"Fe",   26},
    {"Co",   27},
    {"Ni",   28},
    {"Cu",   29},
    {"Zn",   30},
    {"Ga",   31},
    {"Ge",   32},
    {"As",   33},
    {"Se",   34},
    {"Br",   35},
    {"Kr",   36},
    {"Rb",   37},
    {"Sr",   38},
    {"Y",    39},
    {"Zr",   40},
    {"Nb",   41},
    {"Mo",   42},
    {"Tc",   43},
    {"Ru",   44},
    {"Rh",   45},
    {"Pd",   46},
    {"Ag",   47},
    {"Cd",   48},
    {"In",   49},
    {"Sn",   50},
    {"Sb",   51},
    {"Te",   52},
    {"I",    53},
    {"Xe",   54},
    {"Cs",   55},
    {"Ba",   56},
    {"La",   57},
    {"Ce",   58},
    {"Pr",   59},
    {"Nd",   60},
    {"Pm",   61},
    {"Sm",   62},
    {"Eu",   63},
    {"Gd",   64},
    {"Tb",   65},
    {"Dy",   66},
    {"Ho",   67},
    {"Er",   68},
    {"Tm",   69},
    {"Yb",   70},
    {"Lu",   71},
    {"Hf",   72},
    {"Ta",   73},
    {"W",    74},
    {"Re",   75},
    {"Os",   76},
    {"Ir",   77},
    {"Pt",   78},
    {"Au",   79},
    {"Hg",   80},
    {"Tl",   81},
    {"Pb",   82},
    {"Bi",   83},
    {"Po",   84},
    {"At",   85},
    {"Rn",   86},
    {"Fr",   87},
    {"Ra",   88},
    {"Ac",   89},
    {"Th",   90},
    {"Pa",   91},
    {"U",    92},
    {"Np",   93},
    {"Pu",   94},
    {"Am",   95},
    {"Cm",   96},
    {"Bk",   97},
    {"Cf",   98},
    {"Es",   99},
    {"Fm",  100},
    {"Md",  101},
    {"No",  102},
    {"Lr",  103},
    {"Rf",  104},
    {"Db",  105},
    {"Sg",  106},
    {"Bh",  107},
    {"Hs",  108},
    {"Mt",  109},
    {"Ds",  110},
    {"Rg",  111},
    {"Cn",  112},
    {"Uut", 113},
    {"Fl",  114},
    {"Uup", 115},
    {"Lv",  116},
    {"Uus", 117},
    {"Og",  118}
};
// clang-format on


/**
 *  Use the .begin() and .end() iterators of a std::vector to construct a boost::bimap
 *
 *  Adapted from https://stackoverflow.com/a/20290421/7930415
 */
static const boost::bimap<std::string, size_t> periodic_table {elements_list.begin(), elements_list.end()};


/*
 * FUNCTION IMPLEMENTATIONS
 */

/**
 *  @param symbol       the name of an element
 *
 *  @return the atomic number of the corresponding element
 */
size_t elementToAtomicNumber(const std::string& symbol) {

    return periodic_table.left.at(symbol);
}


/**
 *  @param atomic_number    the atomic number of an element
 *
 *  @return the symbol of the corresponding element
 */
const std::string& atomicNumberToElement(size_t atomic_number) {

    return periodic_table.right.at(atomic_number);
}


}  // namespace elements
}  // namespace GQCP
