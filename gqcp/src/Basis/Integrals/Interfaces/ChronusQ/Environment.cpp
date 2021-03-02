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


#include "Basis/Integrals/Interfaces/ChronusQ/Environment.hpp"

#include "Basis/Integrals/Interfaces/ChronusQ/math.hpp"
#include "Basis/Integrals/Interfaces/LibintInterfacer.hpp"


namespace ChronusQ {


std::vector<std::vector<std::array<int, 3>>> Environment::pop_cart_ang_list() {

    std::vector<std::vector<std::array<int, 3>>> list;

    // generate angular momentum list, can only be called once.
    int k, xx, yy, x, y, z;
    for (k = 0; k <= LIBINT2_MAX_AM; k++) {
        list.emplace_back();  //loop over possible angular momentum
        for (xx = 0; xx < k + 1; xx++) {
            x = k - xx;
            for (yy = 0; yy < xx + 1; yy++) {
                y = xx - yy;
                z = k - x - y;
                list[k].push_back({x, y, z});
            }
        }
    }

    return list;
}
std::vector<std::vector<std::array<int, 3>>> Environment::cart_ang_list = Environment::pop_cart_ang_list();


std::vector<std::vector<double>> Environment::pop_car2sph_matrix() {

    std::vector<std::vector<double>> matrix;

    int l[3];
    int m;
    int carsize;  //cartesian size
    double scalecoeff;

    for (int L = 0; L <= LIBINT2_MAX_AM; L++) {

        carsize = (L + 1) * (L + 2) / 2;
        matrix.emplace_back((2 * L + 1) * carsize, 0.0);
        if (L == 0) {
            matrix[L][0] = 1.0;
        } else if (L == 1) {
            matrix[L][0 * 3 + 0] = 1.0;
            matrix[L][1 * 3 + 1] = 1.0;
            matrix[L][2 * 3 + 2] = 1.0;
        } else {
            for (int p = 0; p < L + 1; p++)
                for (int q = 0; q < carsize; q++) {
                    for (int k = 0; k < 3; k++)
                        l[k] = Environment::cart_ang_list[L][q][k];
                    m = -L + p;
                    if (m < 0) {
                        auto cplxcoeff = cart2sphCoeff(L, m, l[0], l[1], l[2]);  // complex coefficient

                        matrix[L][p * carsize + q] = sqrt(2.0) * (-cplxcoeff.imag());

                        matrix[L][(2 * L - p) * carsize + q] = sqrt(2.0) * cplxcoeff.real();

                        double scalecoeff = sqrt(doubleFact(L) /
                                                 (doubleFact(l[0]) * doubleFact(l[1]) * doubleFact(l[2])));
                        matrix[L][p * carsize + q] = matrix[L][p * carsize + q] * scalecoeff;
                        matrix[L][(2 * L - p) * carsize + q] = matrix[L][(2 * L - p) * carsize + q] * scalecoeff;
                    } else if (m == 0) {
                        matrix[L][p * carsize + q] = cart2sphCoeff(L, m, l[0], l[1], l[2]).real();
                        scalecoeff = sqrt(doubleFact(L) /
                                          (doubleFact(l[0]) * doubleFact(l[1]) * doubleFact(l[2])));
                        matrix[L][p * carsize + q] = matrix[L][p * carsize + q] * scalecoeff;
                    }
                }
        }

    }  //loop over angular momentum


    return matrix;
}  //pop_car2sph_matrix()
std::vector<std::vector<double>> Environment::car2sph_matrix = Environment::pop_car2sph_matrix();

std::array<std::array<double, 25>, 3201> Environment::generateFmTTable() {

    std::array<std::array<double, 25>, 3201> table;

    double intervalFmT = 0.025;
    double T = 0.0;
    int MaxTotalL = 25;
    int MaxFmTPt = 3201;
    double critT = 33.0;  // critical value for T. for T>critT, use limit formula
    double expT, factor, term, sum, twoT, Tn;
    for (int i = 0; i < MaxFmTPt; i++) {
        if (std::abs(T) <= 1.0e-10) {
            for (int m = 0; m <= MaxTotalL; m++)
                table[i][m] = 1.0 / (2.0 * m + 1);
        } else if (T > critT) {
            table[i][0] = 0.5 * sqrt(M_PI / T);
            twoT = 2.0 * T;
            Tn = 1.0;
            for (int m = 1; m < MaxTotalL; m++) {
                Tn *= twoT;
                table[i][m] = table[i][m - 1] * (2 * m - 1) / twoT;
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
            table[i][MaxTotalL] = expT * sum;
            twoT = 2.0 * T;
            for (int m = MaxTotalL - 1; m >= 0; m--)
                table[i][m] = (twoT * table[i][m + 1] + expT) / (2 * m + 1);
        }  // else
        T += intervalFmT;
    }  // for i
}  // generateFmTTable()
std::array<std::array<double, 25>, 3201> Environment::FmTTable = Environment::generateFmTTable();


}  // namespace ChronusQ
