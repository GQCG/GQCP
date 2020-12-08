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

#pragma once


#include "Operator/SecondQuantized/ModelHamiltonian/HoppingMatrix.hpp"
#include "Operator/SecondQuantized/RSQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/RSQTwoElectronOperator.hpp"


namespace GQCP {


/**
 *  The Hubbard model Hamiltonian.
 * 
 *  @tparam _Scalar         The scalar type for a hopping matrix element.
 */
template <typename _Scalar>
class HubbardHamiltonian {
public:
    // The scalar type for a hopping matrix element.
    using Scalar = _Scalar;


private:
    // The Hubbard hopping matrix.
    HoppingMatrix<Scalar> H;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  Create a `HubbardHamiltonian` from a `HoppingMatrix`.
     * 
     *  @param H            The Hubbard hopping matrix.
     */
    HubbardHamiltonian(const HoppingMatrix<Scalar>& H) :
        H {H} {}


    /*
     *  MARK: Integral access
     */

    /**
     *  @return The core Hamiltonian (i.e. resulting from the Hubbard hopping operator) as a one-electron operator.
     * 
     *  @note This method is only available for real scalars.
     */
    template <typename Z = Scalar>
    enable_if_t<std::is_same<Z, double>::value, ScalarRSQOneElectronOperator<double>> core() const {

        // Prepare some variables.
        const auto K = this->numberOfLatticeSites();
        const auto& H = this->hoppingMatrix();
        auto h = ScalarRSQOneElectronOperator<double>::Zero(K);

        // The one-electron hopping terms can be found on the off-diagonal elements of the hopping matrix.
        for (size_t p = 0; p < K; p++) {
            for (size_t q = p; q < K; q++) {
                if (p != q) {
                    h.parameters()(p, q) = H(p, q);
                    h.parameters()(q, p) = H(q, p);
                }
            }
        }

        return h;
    }


    /**
     *  @return The two-electron part of the Hamiltonian (resulting from the on-site repulsion) as a two-electron operator.
     * 
     *  @note This method is only available for real scalars.
     */
    template <typename Z = Scalar>
    enable_if_t<std::is_same<Z, double>::value, ScalarRSQTwoElectronOperator<double>> twoElectron() const {

        const auto K = this->numberOfLatticeSites();
        SquareRankFourTensor<double> g_par = SquareRankFourTensor<double>::Zero(K);

        // The two-electron on-site repulsion is found on the diagonal of the hopping matrix.
        for (size_t p = 0; p < K; p++) {
            g_par(p, p, p, p) = this->hoppingMatrix()(p, p);
        }

        return ScalarRSQTwoElectronOperator<double> {g_par};
    }


    /**
     *  @return The Hubbard hopping matrix for this Hubbard model Hamiltonian.
     */
    const HoppingMatrix<Scalar>& hoppingMatrix() const { return this->H; }


    /*
     *  MARK: General information
     */

    /**
     *  @return the number of lattice sites corresponding used in this Hubbard model Hamiltonian
     */
    size_t numberOfLatticeSites() const { return this->H.numberOfLatticeSites(); }
};


}  // namespace GQCP
