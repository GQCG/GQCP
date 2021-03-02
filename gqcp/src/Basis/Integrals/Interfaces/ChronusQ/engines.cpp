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

#include "Basis/Integrals/Interfaces/ChronusQ/engines.hpp"

#include "Basis/Integrals/Interfaces/ChronusQ/Environment.hpp"
#include "Basis/Integrals/Interfaces/ChronusQ/math.hpp"


namespace ChronusQ {


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
