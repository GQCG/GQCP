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


#include "Mathematical/Representation/SquareMatrix.hpp"
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
    // The scalar type for a matrix element.
    using Scalar = _Scalar;


private:
    // The Hubbard Hamiltonian matrix.
    HoppingMatrix<Scalar> hopping_matrix;
    SquareMatrix<Scalar> U_matrix;
    SquareMatrix<Scalar> mu_matrix;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  Create a Hubbard Hamiltonian matrix from its representation as a `SquareMatrix`.
     *
     *  @param H        The Hopping matrix.
     *  @param U        The U contributions, represented as a `SquareMatrix`.
     *  @param mu       The on-site potentials, represented as a `SquareMatrix`.
     */
    HubbardHamiltonian(const HoppingMatrix<Scalar>& H, const SquareMatrix<Scalar>& U, const SquareMatrix<Scalar>& mu) :
        hopping_matrix {H},
        U_matrix {U},
        mu_matrix {mu} {

        const SquareMatrix<Scalar> total = this->hopping_matrix.matrix() + this->U_matrix + this->mu_matrix;

        if (!total.isHermitian()) {
            throw std::invalid_argument("HubbardHamiltonian::HubbardHamiltonian(const HoppingMatrix<Scalar>&, const SquareMatrix<Scalar>&, const SquareMatrix<Scalar>&): The total Hubbard Hamiltonian matrix must be Hermitian.");
        }
    }


    /**
     *  Create a `HubbardHamiltonian` from a `HoppingMatrix` and constant parameters `U`and `mu`.
     *
     *  @param H            The Hubbard hopping matrix.
     *  @param U            The on site repulsion value.
     *  @param mu           The on site potential. Default is zero.
     */
    HubbardHamiltonian(const HoppingMatrix<Scalar>& H, const double& U, const double& mu = 0.0) :
        hopping_matrix {H},
        U_matrix {SquareMatrix<double>::Identity(H.matrix().dimension())},
        mu_matrix {SquareMatrix<double>::Identity(H.matrix().dimension())} {

        // Fill in the on site repulsion on the diagonal.
        for (size_t i = 0; i < H.matrix().dimension(); i++) {
            this->U_matrix(i, i) *= U;
            this->mu_matrix(i, i) *= mu;
        }
    }


    /**
     *  Create a `HubbardHamiltonian` from a `HoppingMatrix` and parameters `U`and `mu` for each site as a vector.
     *
     *  @param H            The Hubbard hopping matrix.
     *  @param U            The on site repulsion values as a vector.
     *  @param mu           The on site potential values as a vector.
     */
    HubbardHamiltonian(const HoppingMatrix<Scalar>& H, const std::vector<double>& U, const std::vector<double>& mu) :
        hopping_matrix {H},
        U_matrix {SquareMatrix<double>::Identity(H.matrix().dimension())},
        mu_matrix {SquareMatrix<double>::Identity(H.matrix().dimension())} {

        // Fill in the given vectors.
        for (size_t i = 0; i < U.size(); i++) {
            this->U_matrix(i, i) *= U[i];
            this->mu_matrix(i, i) *= mu[i];
        }
    }


    /*
     * MARK: Named constructors
     */

    /**
     *  Create a random Hubbard Hamiltonian matrix matrix with elements distributed uniformly in [-1.0, 1.0].
     *
     *  @param K        The number of lattice sites.
     *
     *  @return A random Hubbard Hamiltonian matrix.
     *
     *  @note This method is only available for real scalars.
     *  @note This method always sets the on-site potential contributions (mu) to zero.
     */
    template <typename Z = Scalar>
    static enable_if_t<std::is_same<Z, double>::value, HubbardHamiltonian<double>> Random(const size_t K) {
        // prepare the empty matrices.
        auto mu = SquareMatrix<double>::Zero(K);
        auto hopping = SquareMatrix<double>::Zero(K);
        auto U = SquareMatrix<double>::Zero(K);

        // Generate a random matrix.
        const auto hopping_plus_u = SquareMatrix<double>::RandomSymmetric(K);

        // Distribute the random matrix values correctly between the component matrices.
        for (size_t p = 0; p < K; p++) {
            for (size_t q = p; q < K; q++) {
                if (p != q) {
                    hopping(p, q) = hopping_plus_u(p, q);
                    hopping(q, p) = hopping_plus_u(q, p);
                } else if (p == q) {
                    U(p, q) = hopping_plus_u(p, q);
                }
            }
        }

        return HubbardHamiltonian {HoppingMatrix<double> {hopping}, U, mu};
    }


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
        const auto& H = this->oneElectronContributions();
        auto h = ScalarRSQOneElectronOperator<double>::Zero(K);

        // The one-electron hopping terms can be found on the off-diagonal elements of the hopping matrix.
        for (size_t p = 0; p < K; p++) {
            for (size_t q = p; q < K; q++) {
                h.parameters()(p, q) = H(p, q);
                h.parameters()(q, p) = H(q, p);
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
        const auto& H = this->onSiteRepulsionMatrix();
        SquareRankFourTensor<double> g_par = SquareRankFourTensor<double>::Zero(K);

        // The two-electron on-site repulsion is found on the diagonal of the hopping matrix.
        for (size_t p = 0; p < K; p++) {
            g_par(p, p, p, p) = H(p, p);
        }

        return ScalarRSQTwoElectronOperator<double> {g_par};
    }


    /*
     *  MARK: General information
     */

    /**
     *  @return The Hubbard hopping matrix for this Hubbard model Hamiltonian.
     */
    const HoppingMatrix<Scalar>& hoppingMatrix() const { return this->hopping_matrix; }


    /**
     *  @return the number of lattice sites corresponding used in this Hubbard model Hamiltonian
     */
    size_t numberOfLatticeSites() const { return this->hopping_matrix.matrix().dimension(); }


    /**
     *  @return The one electron parameter contributions of this Hubbard model Hamiltonian.
     */
    const SquareMatrix<Scalar> oneElectronContributions() const {
        return this->hopping_matrix.matrix() + this->mu_matrix;
    }


    /**
     *  @return The matrix containing the on-site potentials for this Hubbard model Hamiltonian.
     */
    const SquareMatrix<Scalar>& onSitePotentialMatrix() const { return this->mu_matrix; }


    /**
     *  @return The matrix containing the on-site repulsions for this Hubbard model Hamiltonian.
     */
    const SquareMatrix<Scalar>& onSiteRepulsionMatrix() const { return this->U_matrix; }
};


}  // namespace GQCP
