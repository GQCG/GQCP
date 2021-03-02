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


#include "Basis/Integrals/Interfaces/ChronusQ/math.hpp"

#include "Basis/Integrals/Interfaces/ChronusQ/Environment.hpp"

#include <iostream>


namespace ChronusQ {


double factorial(int t) {
    int i;
    double tmp = 1.0;
    if (t < 0)
        std::cout << "Factorial (t!) only defined on domain t in [0,inf)" << std::endl;
    if (t == 0)
        return 1.0;
    else {
        for (i = 1; i <= t; i++)
            tmp *= i;
        return tmp;
    }
}  //factorial


double doubleFact(int t) {
    int i;
    double tmp = 1.0;
    if (t < 0)
        std::cout << "Double factorial (t!!) only defined on domain t in [0,inf)"
                  << std::endl;
    if (t == 0)
        return 1.0;
    else {
        for (i = 1; i <= t; i++)
            tmp *= (2 * i - 1);
        return tmp;
    }
}  // doubleFact


double polyCoeff(int l, int i) {
    if (l >= i)
        return factorial(l) / (factorial(i) * factorial(l - i));
    else
        std::cout << "polyCoeff error" << std::endl;
}


std::complex<double> cart2sphCoeff(int L, int m, int lx, int ly, int lz) {

    //calculate the cartesian to spherical transformation coefficient.

    int Ltotal;
    dcomplex coeff(0.0);
    Ltotal = lx + ly + lz;
    double tmp = 0.0;
    if (L != Ltotal) {
        return coeff;
    }
    double j;
    j = (double(lx + ly) - std::abs(double(m))) / 2;
    if (fmod(j, 1) > 0) {
        return coeff;
    }
    dcomplex sumval(0.0);
    dcomplex ttmmpp, sumsumval;
    dcomplex pref, absmchooselxm2k, ichoosej;
    int i, k;
    if (Ltotal == L) {
        pref = sqrt(factorial(lx * 2) * factorial(2 * ly) * factorial(2 * lz) * factorial(L) * factorial(L - std::abs(m)) / (factorial(2 * L) * factorial(lx) * factorial(ly) * factorial(lz) * factorial(L + std::abs(m)))) / (factorial(L) * pow(2, L));

        i = 0;

        while (i <= double((L - std::abs(m)) / 2)) {
            sumsumval = 0.0;
            for (k = 0; k <= j; k++) {
                if (m >= 0) {
                    ttmmpp = double(std::abs(m) - lx + 2 * k) / 2;
                } else {
                    ttmmpp = -double(std::abs(m) - lx + 2 * k) / 2;
                }

                if ((std::abs(m) >= (lx - 2 * k)) && ((lx - 2 * k) >= 0)) {
                    absmchooselxm2k = polyCoeff(std::abs(m), lx - 2 * k);
                } else {
                    absmchooselxm2k = 0.0;
                }
                sumsumval = sumsumval + polyCoeff(j, k) * absmchooselxm2k * pow(-1.0, ttmmpp);
            }
            if (i < j || (j < 0)) {
                ichoosej = 0.0;
            } else {
                ichoosej = polyCoeff(i, j);
            }
            sumval = sumval + polyCoeff(L, i) * ichoosej * pow(-1, i) * factorial(2 * L - 2 * i) /
                                  (factorial(L - std::abs(m) - 2 * i)) * sumsumval;
            i = i + 1;
        }
        coeff = pref * sumval;
        return coeff;
    }
}  // cart2sphCoeff


void cart2sph_complex_transform(int l_i, int l_j, std::vector<dcomplex>& shell_element_sph, std::vector<dcomplex>& shell_element_cart) {

    int cart_i = (l_i + 1) * (l_i + 2) / 2;
    int cart_j = (l_j + 1) * (l_j + 2) / 2;
    int cartsize = cart_i * cart_j;
    int sphsize = (2 * l_i + 1) * (2 * l_j + 1);
    dcomplex tempVal;

    if (sphsize != shell_element_sph.size())
        std::cout << "spherical dimension doesn't match" << std::endl;
    if (cartsize != shell_element_cart.size())
        std::cout << "cartesian dimension doesn't match" << std::endl;

    for (int i = 0; i < 2 * l_i + 1; i++) {
        for (int j = 0; j < 2 * l_j + 1; j++) {
            tempVal = 0.0;
            for (int p = 0; p < cart_i; p++) {
                for (int q = 0; q < cart_j; q++) {
                    tempVal += Environment::car2sph_matrix[l_i][i * cart_i + p] * Environment::car2sph_matrix[l_j][j * cart_j + q] * shell_element_cart[p * cart_j + q];
                }
            }

            if (std::abs(tempVal) < 1.0e-15)
                tempVal = 0.0;
            shell_element_sph[i * (2 * l_j + 1) + j] = tempVal;
        }
    }
}  //cart2sph_complex_transform


}  // namespace ChronusQ
