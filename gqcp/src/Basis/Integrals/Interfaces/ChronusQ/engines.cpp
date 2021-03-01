/*
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *
 *  Copyright (C) 2014-2018 Li Research Group (University of Washington)
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
 *
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

#include "Basis/Integrals/Interfaces/ChronusQ/engines.hpp"


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


std::vector<std::vector<dcomplex>> ComplexGIAOIntEngine::computeGIAOOverlapS(
    libint2::ShellPair& pair, libint2::Shell& shell1, libint2::Shell& shell2, const std::array<double, 3>& H) {


    int nElement = Environment::cart_ang_list[shell1.contr[0].l].size() * Environment::cart_ang_list[shell2.contr[0].l].size();

    // evaluate the number of integral in the shell pair with cartesian gaussian
    std::vector<dcomplex> S_cartshell;

    dcomplex S;
    int lA[3], lB[3];

    double ka[3], kb[3];
    // here we define k as the direct phase factor without negative sign:
    // w = G* exp[ik dot (r-R)] for both bra and ket
    // ka (bra) = -1/2  A X H
    // kb (ket) = 1/2 B X H

    /*
    for ( int i = 0 ; i < 3 ; i++ ){
      ka[i] = 
      kb[i] = 
    } 
*/
    // std::cout<<"H "<<H[0]<<" "<<H[1]<<" "<<H[2]<<std::endl;

    ka[0] = -0.5 * (shell1.O[1] * H[2] - shell1.O[2] * H[1]);
    ka[1] = -0.5 * (shell1.O[2] * H[0] - shell1.O[0] * H[2]);
    ka[2] = -0.5 * (shell1.O[0] * H[1] - shell1.O[1] * H[0]);

    kb[0] = 0.5 * (shell2.O[1] * H[2] - shell2.O[2] * H[1]);
    kb[1] = 0.5 * (shell2.O[2] * H[0] - shell2.O[0] * H[2]);
    kb[2] = 0.5 * (shell2.O[0] * H[1] - shell2.O[1] * H[0]);


    //std::cout<<"ka "<<ka[0]<<"  "<<ka[1]<<"  "<<ka[2]<<std::endl;
    //std::cout<<"kb "<<kb[0]<<"  "<<kb[1]<<"  "<<kb[2]<<std::endl;

    double K[3];
    for (int mu = 0; mu < 3; mu++)
        K[mu] = ka[mu] + kb[mu];

    auto ss_shellpair = computecompOverlapss(pair, shell1, ka, shell2, kb);


    //for ( auto sselement : ss_shellpair ){
    //  std::cout<<"SS in a shellpair "<<std::setprecision(12)<<sselement<<std::endl;
    //}

    for (int i = 0; i < Environment::cart_ang_list[shell1.contr[0].l].size(); i++) {

        for (int j = 0; j < Environment::cart_ang_list[shell2.contr[0].l].size(); j++) {
            for (int k = 0; k < 3; k++) {

                lA[k] = Environment::cart_ang_list[shell1.contr[0].l][i][k];
                lB[k] = Environment::cart_ang_list[shell2.contr[0].l][j][k];
            }


            S = comphRRSab(pair, shell1, shell2, K, ss_shellpair,
                           shell1.contr[0].l, lA, shell2.contr[0].l, lB);

            // std::cout<<"S value"<<std::setprecision(12)<<S<<std::endl;

            S_cartshell.push_back(S);

        }  // loop over ij
    }

    if ((not shell1.contr[0].pure) and (not shell2.contr[0].pure)) {
        // if both sides are cartesian, return cartesian gaussian integrals

        return {S_cartshell};
    }

    std::vector<std::vector<dcomplex>> S_shellpair_sph(1);

    S_shellpair_sph[0].assign(((2 * shell1.contr[0].l + 1) * (2 * shell2.contr[0].l + 1)), 0.0);

    cart2sph_complex_transform(shell1.contr[0].l, shell2.contr[0].l,
                               S_shellpair_sph[0], S_cartshell);

    return S_shellpair_sph;


}  // computeGIAOOverlapS


// compute uncontracted overlap of (s||s) type for a shellpair

std::vector<dcomplex> ComplexGIAOIntEngine::computecompOverlapss(
    libint2::ShellPair& pair, libint2::Shell& shell1, double* ka, libint2::Shell& shell2, double* kb) {


    std::vector<dcomplex> ss_shellpair;
    dcomplex onei;
    onei.real(0);
    onei.imag(1);
    for (auto& pripair : pair.primpairs) {

        dcomplex norm;
        dcomplex tmpVal;
        norm = shell1.contr[0].coeff[pripair.p1] * shell2.contr[0].coeff[pripair.p2];

        double realpart = 0.0;
        for (int mu = 0; mu < 3; mu++) {
            realpart -= pow((ka[mu] + kb[mu]), 2);
        }  // for mu

        realpart *= 0.25 * pripair.one_over_gamma;

        double imagpart = 0.0;
        for (int mu = 0; mu < 3; mu++) {
            imagpart += ka[mu] * (pripair.P[mu] - shell1.O[mu]) + kb[mu] * (pripair.P[mu] - shell2.O[mu]);
        }  // for mu

        dcomplex z = realpart + imagpart * onei;

        tmpVal = norm * exp(z) *
                 pow(sqrt(M_PI), 3) * sqrt(pripair.one_over_gamma) * pripair.K;

        ss_shellpair.push_back(tmpVal);
    }

    return ss_shellpair;
}


/**
   *  \brief Perform the horizontal recurrence relation for the contracted overlap integral
   *
   *  (a|b) = (A-B)(a|b-1) + (a+1|b-1)
   *
   *  where a,b are the angular momentum, A,B are the nuclear coordinates.
   *
   *  \param [in] pair    Shell pair data for shell1, shell2
   *  \param [in] shell1  Bra shell
   *  \param [in] shell2  Ket shell
   *  \param [in] LA      total Bra angular momentum
   *  \param [in] lA      Bra angular momentum vector (lAx,lAy,lAz)
   *  \param [in] LB      total Ket angular momentum
   *  \param [in] lB      Ket angular momentum vector (lBx,lBy,lBz)
   *
   *  \returns a contracted overlap integral
   *
   */
dcomplex ComplexGIAOIntEngine::comphRRSab(libint2::ShellPair& pair, libint2::Shell& shell1,
                                          libint2::Shell& shell2, double* K, std::vector<dcomplex>& ss_shellpair,
                                          int LA, int* lA, int LB, int* lB) {


    int iWork, lAp1[3], lBm1[3];
    dcomplex tmpVal = 0.0;

    if (LB > LA) {  // if LB>LA, use horizontal recursion to make LA>LB.

        for (iWork = 0; iWork < 3; iWork++) {
            lAp1[iWork] = lA[iWork];
            lBm1[iWork] = lB[iWork];
        };

        if (lB[0] > 0)
            iWork = 0;
        else if (lB[1] > 0)
            iWork = 1;
        else if (lB[2] > 0)
            iWork = 2;
        lAp1[iWork]++;
        lBm1[iWork]--;

        tmpVal = comphRRSab(pair, shell1, shell2, K, ss_shellpair,
                            LA + 1, lAp1, LB - 1, lBm1);
        tmpVal += pair.AB[iWork] * comphRRSab(pair, shell1, shell2, K,
                                              ss_shellpair, LA, lA, LB - 1, lBm1);

        return tmpVal;

    }  // if ( LB > LA )

    if (LA == 0) {
        // (s|s)
        auto pripairindex = 0;
        for (auto& pripair : pair.primpairs) {
            tmpVal += ss_shellpair[pripairindex];
            pripairindex++;

            /*
        double norm;
        norm = shell1.contr[0].coeff[pripair.p1]* shell2.contr[0].coeff[pripair.p2];
    
        double realpart=0.0;
        for ( int mu = 0 ; mu < 3 ; mu++ ) {
          realpart -= pow( (ka[mu]+kb[mu]), 2 );
        } // for mu 
        realpart *= 0.25*pripair.one_over_gamma; 
 
        double imagpart=0.0;
        for ( int mu = 0 ; mu < 3 ; mu++ ) {
          imagpart += ka[mu]*(pripair.P[mu] - shell1.O[mu]) 
                      + kb[mu]*(pripair.P[mu] - shell2.O[mu]);
        }  // for mu 

        dcomplex z = ( realpart , imagpart );   
        
        tmpVal += norm * exp ( z ) *
                  pow(sqrt(M_PI),3) * sqrt(pripair.one_over_gamma)*pripair.K ;
*/

        }  // for pripair

        return tmpVal;
    } else if (LB == 0) {
        // (|s)
        auto pripairindex = 0;
        for (auto pripair : pair.primpairs) {

            //tmpVal+= shell1.contr[0].coeff[pripair.p1]* shell2.contr[0].coeff[pripair.p2]*


            tmpVal +=
                compvRRSa0(pripair, shell1, K, ss_shellpair[pripairindex], LA, lA);

            pripairindex++;

        }  // for pripair

        return tmpVal;
    };  // else if(LB == 0)

    // here LB > 0

    for (iWork = 0; iWork < 3; iWork++) {
        lAp1[iWork] = lA[iWork];
        lBm1[iWork] = lB[iWork];
    };
    if (lB[0] > 0)
        iWork = 0;
    else if (lB[1] > 0)
        iWork = 1;
    else if (lB[2] > 0)
        iWork = 2;
    lAp1[iWork]++;
    lBm1[iWork]--;


    tmpVal = comphRRSab(pair, shell1, shell2, K, ss_shellpair, LA + 1, lAp1, LB - 1, lBm1);
    tmpVal += (shell1.O[iWork] - shell2.O[iWork]) * comphRRSab(pair, shell1, shell2, K, ss_shellpair, LA, lA, LB - 1, lBm1);

    //  tmpVal+= pair.AB[iWork]*hRRSab(pair,shell1,shell2,LA,lA,LB-1,lBm1);
    return tmpVal;

}  // comphRRSab


//
//------------------------------------//
// overlap horizontal recursion iPP   //
// (a|b) = (A-B)(a|b-1) + (a+1|b-1)   //
//------------------------------------//

/**
   *  \brief Perform the horizontal recurrence relation for the contracted overlap integral
   *
   *  (a|b) = (A-B)(a|b-1) + (a+1|b-1)
   *
   *  where a,b are the angular momentum, A,B are the nuclear coordinates.
   *
   *  \param [in] pair    Shell pair data for shell1, shell2
   *  \param [in] shell1  Bra shell
   *  \param [in] shell2  Ket shell
   *  \param [in] LA      total Bra angular momentum
   *  \param [in] lA      Bra angular momentum vector (lAx,lAy,lAz)
   *  \param [in] LB      total Ket angular momentum
   *  \param [in] lB      Ket angular momentum vector (lBx,lBy,lBz)
   *
   *  \returns a contracted overlap integral
   *
   */
dcomplex ComplexGIAOIntEngine::comphRRiPPSab(libint2::ShellPair::PrimPairData& pripair,
                                             libint2::Shell& shell1, libint2::Shell& shell2, double* K,
                                             dcomplex sspri, int LA, int* lA, int LB, int* lB) {

    int iWork, lAp1[3], lBm1[3];
    dcomplex tmpVal = 0.0;

    if (LB > LA) {  // if LB>LA, use horizontal recursion to make LA>=LB.

        for (int mu = 0; mu < 3; mu++) {
            lAp1[mu] = lA[mu];
            lBm1[mu] = lB[mu];
        }  // for mu

        if (lB[0] > 0)
            iWork = 0;
        else if (lB[1] > 0)
            iWork = 1;
        else if (lB[2] > 0)
            iWork = 2;
        lAp1[iWork]++;
        lBm1[iWork]--;

        tmpVal = comphRRiPPSab(pripair, shell1, shell2, K, sspri,
                               LA + 1, lAp1, LB - 1, lBm1);
        tmpVal += (shell1.O[iWork] - shell2.O[iWork]) * comphRRiPPSab(pripair, shell1, shell2, K,
                                                                      sspri, LA, lA, LB - 1, lBm1);

        return tmpVal;

    }  // if ( LB > LA )

    if (LA == 0) {
        // (s||s)
        tmpVal = sspri;

        return tmpVal;

    } else if (LB == 0) {
        // (|s)

        tmpVal += compvRRSa0(pripair, shell1, K, sspri, LA, lA);

        return tmpVal;
    };  // else if(LB == 0)

    // here LB > 0

    for (iWork = 0; iWork < 3; iWork++) {
        lAp1[iWork] = lA[iWork];
        lBm1[iWork] = lB[iWork];
    }  // for iWork

    if (lB[0] > 0)
        iWork = 0;
    else if (lB[1] > 0)
        iWork = 1;
    else if (lB[2] > 0)
        iWork = 2;
    lAp1[iWork]++;
    lBm1[iWork]--;

    tmpVal = comphRRiPPSab(pripair, shell1, shell2, K, sspri, LA + 1, lAp1, LB - 1, lBm1);
    tmpVal += (shell1.O[iWork] - shell2.O[iWork]) * comphRRiPPSab(pripair, shell1, shell2, K, sspri, LA, lA, LB - 1, lBm1);
    //  tmpVal+= pair.AB[iWork]*hRRSab(pair,shell1,shell2,LA,lA,LB-1,lBm1);
    return tmpVal;

}  // comphRRiPPSab


//----------------------------------------------------------//
// complex overlap vertical recursion                               //
// (a|0) = (P-A)(a-1|0) + halfInvZeta*N_(a-1)*(a-2|0)       //
//----------------------------------------------------------//
/**
   *  \brief Perform the vertical recurrence relation for the uncontracted overlap integral 
   *
   *  (a|0) = (P-A)(a-1|0) + 1/2 *1/Zeta * N_(a-1)*(a-2|0)
   *
   *  where a is angular momentum, Zeta=zeta_a+zeta_b, A is bra nuclear coordinate. 
   *  P = (zeta_a*A+zeta_b*B)/Zeta
   *
   *  \param [in] pripair Primitive Shell pair data for shell1, shell2
   *  \param [in] shell1  Bra shell
   *  \param [in] LA      total Bra angular momentum
   *  \param [in] lA      Bra angular momentum vector (lAx,lAy,lAz)
   *
   *  \returns an uncontracted overlap integral
   *
   */
dcomplex ComplexGIAOIntEngine::compvRRSa0(libint2::ShellPair::PrimPairData& pripair,
                                          libint2::Shell& shell1, double* K, dcomplex sspri, int LA, int* lA) {

    //notice: contraction coeffs are included in sspri.
    dcomplex tmpVal = 0.0;

    if (LA == 0) {  //[s||s]

        tmpVal = sspri;
        // tmpVal = pow(sqrt(M_PI),3) * sqrt(pripair.one_over_gamma)*pripair.K ;
        return tmpVal;
    }  // if (LA == 0)

    int iWork;
    int lAm1[3];
    for (iWork = 0; iWork < 3; iWork++)
        lAm1[iWork] = lA[iWork];
    if (lA[0] > 0)
        iWork = 0;
    else if (lA[1] > 0)
        iWork = 1;
    else if (lA[2] > 0)
        iWork = 2;

    /*
    if(LA == 1) {
     tmpVal = (pripair.P[iWork]-shell1.O[iWork])*pow(sqrt(M_PI),3) * sqrt(pripair.one_over_gamma)*pripair.K ; 
     return tmpVal;
    }
  */

    lAm1[iWork]--;
    dcomplex onei;
    onei.real(0);
    onei.imag(1);
    tmpVal = (pripair.P[iWork] - shell1.O[iWork] + 0.5 * onei * pripair.one_over_gamma * K[iWork]) * compvRRSa0(pripair, shell1, K, sspri, LA - 1, lAm1);

    //  if(LA == 2&&lA[iWork] == 2)
    //    tmpVal += 1/2*pripair.one_over_gamma
    //            *pow(sqrt(M_PI),3) * sqrt(pripair.one_over_gamma)*pripair.K;
    //  else if (lA[iWork] >=2) {
    if (lA[iWork] >= 2) {
        lAm1[iWork]--;
        tmpVal += (lA[iWork] - 1) * 0.5 * pripair.one_over_gamma * compvRRSa0(pripair, shell1, K, sspri, LA - 2, lAm1);
        //  }
    }  // if ( lA[iWork] >=2 )
    return tmpVal;

}  // compvRRSa0


}  // namespace ChronusQ
