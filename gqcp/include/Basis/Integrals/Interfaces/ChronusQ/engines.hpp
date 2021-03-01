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


// FIXME: Licencing

#pragma once

#include "Basis/Integrals/Interfaces/LibintInterfacer.hpp"

#include <vector>


namespace ChronusQ {

using dcomplex = std::complex<double>;


struct Environment {

    static std::vector<std::vector<std::array<int, 3>>> cart_ang_list;
    static std::vector<std::vector<std::array<int, 3>>> pop_cart_ang_list();

    static std::vector<std::vector<double>> car2sph_matrix;
    static std::vector<std::vector<double>> pop_car2sph_matrix();


    static std::array<std::array<double, 25>, 3201> FmTTable;
    static std::array<std::array<double, 25>, 3201> generateFmTTable();
};


double factorial(int);
double doubleFact(int);
double polyCoeff(int, int);


std::complex<double> cart2sphCoeff(int, int, int, int, int);
void cart2sph_complex_transform(int, int, std::vector<dcomplex>&, std::vector<dcomplex>&);


struct ComplexGIAOIntEngine {

    static std::vector<std::vector<dcomplex>> computeGIAOOverlapS(libint2::ShellPair&, libint2::Shell&, libint2::Shell&, const std::array<double, 3>&);

    // calculate the uncontracted overlap of (s||s) type for a shellpair
    static std::vector<dcomplex> computecompOverlapss(libint2::ShellPair&, libint2::Shell&, double*, libint2::Shell&, double*);

    // complex overlap horizontal recursion for contracted case
    static dcomplex comphRRSab(libint2::ShellPair&, libint2::Shell&, libint2::Shell&, double*, std::vector<dcomplex>&, int, int*, int, int*);

    // complex overlap horizontal recursion iPP specific for uncontracted case
    static dcomplex comphRRiPPSab(libint2::ShellPair::PrimPairData&, libint2::Shell&, libint2::Shell&, double*, dcomplex, int, int*, int, int*);

    // complex overlap vertical recursion for uncontracted case
    static dcomplex compvRRSa0(libint2::ShellPair::PrimPairData&, libint2::Shell&, double*, dcomplex, int, int*);
};


}  // namespace ChronusQ
