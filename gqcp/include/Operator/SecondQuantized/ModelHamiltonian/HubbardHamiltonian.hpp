// This file is part of GQCG-gqcp.
//
// Copyright (C) 2017-2019  the GQCG developers
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
#pragma once


#include "Operator/SecondQuantized/ModelHamiltonian/HoppingMatrix.hpp"
#include "Operator/SecondQuantized/SQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/SQTwoElectronOperator.hpp"


namespace GQCP {


/**
 *  The Hubbard model Hamiltonian.
 */
template <typename _Scalar>
class HubbardHamiltonian {
public:
    using Scalar = _Scalar;


private:
    HoppingMatrix<Scalar> H;  // the Hubbard hopping matrix


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param H            the Hubbard hopping matrix
     */
    HubbardHamiltonian(const HoppingMatrix<Scalar>& H) :
        H {H} {}


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return the core Hamiltonian (i.e. resulting from the Hubbard hopping operator) as a one-electron operator
     * 
     *  @note This method is only available for real scalars.
     */
    template <typename Z = Scalar>
    enable_if_t<std::is_same<Z, double>::value, ScalarSQOneElectronOperator<Scalar>> core() const {

        const auto K = this->numberOfLatticeSites();
        QCMatrix<double> h_par = QCMatrix<double>::Zero(K, K);  // 'par' for 'parameters'

        // The one-electron hopping terms can be found on the off-diagonal elements of the hopping matrix.
        for (size_t i = 0; i < K; i++) {
            for (size_t j = i; j < K; j++) {
                if (i != j) {
                    h_par(i, j) = this->H(i, j);
                    h_par(j, i) = this->H(j, i);
                }
            }
        }

        return ScalarSQOneElectronOperator<double> {h_par};
    }


    /**
     *  @return the Hubbard hopping matrix for this Hubbard model Hamiltonian
     */
    const HoppingMatrix<Scalar>& hoppingMatrix() const { return this->H; }


    /**
     *  @return the number of lattice sites corresponding used in this Hubbard model Hamiltonian
     */
    size_t numberOfLatticeSites() const { return this->H.numberOfLatticeSites(); }


    /**
     *  @return the two-electron part of the Hamiltonian (resulting from the on-site repulsion) as a two-electron operator
     * 
     *  @note This method is only available for real scalars.
     */
    template <typename Z = Scalar>
    enable_if_t<std::is_same<Z, double>::value, ScalarSQTwoElectronOperator<Scalar>> twoElectron() const {

        const auto K = this->numberOfLatticeSites();

        QCRankFourTensor<double> g_par(K);
        g_par.setZero();

        // The two-electron on-site repulsion is found on the diagonal of the hopping matrix.
        for (size_t i = 0; i < K; i++) {
            for (size_t j = i; j < K; j++) {
                if (i == j) {
                    g_par(i, i, i, i) = H(i, i);
                }
            }
        }

        return ScalarSQTwoElectronOperator<double> {g_par};
    }
};


}  // namespace GQCP
