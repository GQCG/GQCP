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

double factorial(int);
double doubleFact(int);
double polyCoeff(int, int);


std::complex<double> cart2sphCoeff(int, int, int, int, int);
void cart2sph_complex_transform(int, int, std::vector<dcomplex>&, std::vector<dcomplex>&);


static std::vector<std::vector<std::array<int, 3>>> cart_ang_list;
void pop_cart_ang_list() {
    // generate angular momentum list, can only be called once.
    int k, xx, yy, x, y, z;
    for (k = 0; k <= LIBINT2_MAX_AM; k++) {
        cart_ang_list.emplace_back();  //loop over possible angular momentum
        for (xx = 0; xx < k + 1; xx++) {
            x = k - xx;
            for (yy = 0; yy < xx + 1; yy++) {
                y = xx - yy;
                z = k - x - y;
                cart_ang_list[k].push_back({x, y, z});
            }
        }
    }
}

static std::vector<std::vector<double>> car2sph_matrix;
void pop_car2sph_matrix() {
    int l[3];
    int m;
    int carsize;  //cartesian size
    double scalecoeff;

    for (int L = 0; L <= LIBINT2_MAX_AM; L++) {

        carsize = (L + 1) * (L + 2) / 2;
        car2sph_matrix.emplace_back((2 * L + 1) * carsize, 0.0);
        if (L == 0) {
            car2sph_matrix[L][0] = 1.0;
        } else if (L == 1) {
            car2sph_matrix[L][0 * 3 + 0] = 1.0;
            car2sph_matrix[L][1 * 3 + 1] = 1.0;
            car2sph_matrix[L][2 * 3 + 2] = 1.0;
        } else {
            for (int p = 0; p < L + 1; p++)
                for (int q = 0; q < carsize; q++) {
                    for (int k = 0; k < 3; k++)
                        l[k] = cart_ang_list[L][q][k];
                    m = -L + p;
                    if (m < 0) {
                        auto cplxcoeff = cart2sphCoeff(L, m, l[0], l[1], l[2]);  // complex coefficient

                        car2sph_matrix[L][p * carsize + q] = sqrt(2.0) * (-cplxcoeff.imag());

                        car2sph_matrix[L][(2 * L - p) * carsize + q] = sqrt(2.0) * cplxcoeff.real();

                        double scalecoeff = sqrt(doubleFact(L) /
                                                 (doubleFact(l[0]) * doubleFact(l[1]) * doubleFact(l[2])));
                        car2sph_matrix[L][p * carsize + q] = car2sph_matrix[L][p * carsize + q] * scalecoeff;
                        car2sph_matrix[L][(2 * L - p) * carsize + q] = car2sph_matrix[L][(2 * L - p) * carsize + q] * scalecoeff;
                    } else if (m == 0) {
                        car2sph_matrix[L][p * carsize + q] = cart2sphCoeff(L, m, l[0], l[1], l[2]).real();
                        scalecoeff = sqrt(doubleFact(L) /
                                          (doubleFact(l[0]) * doubleFact(l[1]) * doubleFact(l[2])));
                        car2sph_matrix[L][p * carsize + q] = car2sph_matrix[L][p * carsize + q] * scalecoeff;
                    }
                }
        }

    }  //loop over angular momentum

}  //pop_car2sph_matrix()

static std::array<std::array<double, 25>, 3201> FmTTable;

void generateFmTTable() {

    double intervalFmT = 0.025;
    double T = 0.0;
    int MaxTotalL = 25;
    int MaxFmTPt = 3201;
    double critT = 33.0;  // critical value for T. for T>critT, use limit formula
    double expT, factor, term, sum, twoT, Tn;
    for (int i = 0; i < MaxFmTPt; i++) {
        if (std::abs(T) <= 1.0e-10) {
            for (int m = 0; m <= MaxTotalL; m++)
                FmTTable[i][m] = 1.0 / (2.0 * m + 1);
        } else if (T > critT) {
            FmTTable[i][0] = 0.5 * sqrt(M_PI / T);
            twoT = 2.0 * T;
            Tn = 1.0;
            for (int m = 1; m < MaxTotalL; m++) {
                Tn *= twoT;
                FmTTable[i][m] = FmTTable[i][m - 1] * (2 * m - 1) / twoT;
            }
        } else {
            expT = exp(-T);
            factor = MaxTotalL + 0.5;
            term = 0.5 / factor;
            sum = term;
            while (term > 1.0e-10) {
                factor += 1.0;
                term *= T / factor;
                sum += term;
            };
            FmTTable[i][MaxTotalL] = expT * sum;
            twoT = 2.0 * T;
            for (int m = MaxTotalL - 1; m >= 0; m--)
                FmTTable[i][m] = (twoT * FmTTable[i][m + 1] + expT) / (2 * m + 1);
        }  // else
        T += intervalFmT;
    }  // for i
}  // generateFmTTable()


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
