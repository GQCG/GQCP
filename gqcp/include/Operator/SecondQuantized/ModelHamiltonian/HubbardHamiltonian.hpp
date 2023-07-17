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
    SquareMatrix<Scalar> Hubbard_Hamiltonian_matrix;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  Create a Hubbard Hamiltonian matrix from its representation as a `SquareMatrix`.
     *
     *  @param H        The Hubbard Hamiltonian matrix, represented as a `SquareMatrix`.
     */
    HubbardHamiltonian(const SquareMatrix<Scalar>& H) :
        Hubbard_Hamiltonian_matrix {H} {

        if (!H.isHermitian()) {
            throw std::invalid_argument("HubbardHamiltonian::HubbardHamiltonian(const SquareMatrix<Scalar>&): The given matrix must be Hermitian.");
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
        Hubbard_Hamiltonian_matrix {SquareMatrix<double>::Identity(H.matrix().dimension())} {
        // Fill in the on site repulsion on the diagonal.
        for (size_t i = 0; i < H.matrix().dimension(); i++) {
            this->Hubbard_Hamiltonian_matrix(i, i) *= U;
            this->Hubbard_Hamiltonian_matrix(i, i) += mu;
        }

        this->Hubbard_Hamiltonian_matrix += H.matrix();
    }


    /**
     *  Create a `HubbardHamiltonian` from a `HoppingMatrix` and parameters `U`and `mu` for each site as a vector.
     *
     *  @param H            The Hubbard hopping matrix.
     *  @param U            The on site repulsion values as a vector.
     *  @param mu           The on site potential values as a vector.
     */
    HubbardHamiltonian(const HoppingMatrix<Scalar>& H, const std::vector<double>& U, const std::vector<double>& mu) :
        Hubbard_Hamiltonian_matrix {SquareMatrix<double>::Identity(H.matrix().dimension())} {
        // Fill in the given vectors.
        for (size_t i = 0; i < U.size(); i++) {
            this->Hubbard_Hamiltonian_matrix(i, i) *= U[i];
            this->Hubbard_Hamiltonian_matrix(i, i) += mu[i];
        }

        this->Hubbard_Hamiltonian_matrix += H.matrix();
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
     */
    template <typename Z = Scalar>
    static enable_if_t<std::is_same<Z, double>::value, HubbardHamiltonian<double>> Random(const size_t K) {
        return HubbardHamiltonian {SquareMatrix<double>::RandomSymmetric(K)};
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
        const auto& H = this->HubbardHamiltonianMatrix();
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
        const auto& H = this->HubbardHamiltonianMatrix();
        SquareRankFourTensor<double> g_par = SquareRankFourTensor<double>::Zero(K);

        // The two-electron on-site repulsion is found on the diagonal of the hopping matrix.
        for (size_t p = 0; p < K; p++) {
            g_par(p, p, p, p) = H(p, p);
        }

        return ScalarRSQTwoElectronOperator<double> {g_par};
    }


    /**
     *  @return The Hubbard hopping matrix for this Hubbard model Hamiltonian.
     */
    const SquareMatrix<Scalar>& HubbardHamiltonianMatrix() const { return this->Hubbard_Hamiltonian_matrix; }


    /*
     *  MARK: General information
     */

    /**
     *  @return the number of lattice sites corresponding used in this Hubbard model Hamiltonian
     */
    size_t numberOfLatticeSites() const { return this->Hubbard_Hamiltonian_matrix.dimension(); }
};


}  // namespace GQCP
