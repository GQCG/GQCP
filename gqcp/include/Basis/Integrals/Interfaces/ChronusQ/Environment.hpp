/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2020 Li Research Group (University of Washington)
 *  
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *  
 *  Contact the Developers:
 *    E-Mail: xsli@uw.edu
 */

// The interfaces and implementations from this file have been adapted from
// those in the Chronus Quantum (ChronusQ) software package.


#pragma once


#include <array>
#include <vector>


namespace ChronusQ {


/**
 *  A wrapper around ChronusQ's global variables.
 */
struct Environment {

    static std::vector<std::vector<std::array<int, 3>>> cart_ang_list;
    static std::vector<std::vector<std::array<int, 3>>> pop_cart_ang_list();

    static std::vector<std::vector<double>> car2sph_matrix;
    static std::vector<std::vector<double>> pop_car2sph_matrix();


    static std::array<std::array<double, 25>, 3201> FmTTable;
    static std::array<std::array<double, 25>, 3201> generateFmTTable();
};


}  // namespace ChronusQ
