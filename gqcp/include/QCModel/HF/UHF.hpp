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


#include "Basis/TransformationMatrix.hpp"
#include "Mathematical/Representation/Matrix.hpp"
#include "QCModel/HF/RHF.hpp"


namespace GQCP {
namespace QCModel {


/**
 *  The unrestricted Hartree-Fock wave function model.
 * 
 *  @tparam _Scalar         the scalar representation of the coefficients in the coefficient matrix
 */
template <typename _Scalar>
class UHF {
public:
    using Scalar = _Scalar;


private:
    size_t N_alpha;  // the number of alpha electrons
    size_t N_beta;   // the number of beta electrons

    VectorX<double> orbital_energies_alpha;  // sorted by ascending energy
    VectorX<double> orbital_energies_beta;   // sorted by ascending energy

    TransformationMatrix<Scalar> C_alpha;  // the coefficient matrix that expresses every alpha spatial orbital (as a column) in its underlying scalar basis
    TransformationMatrix<Scalar> C_beta;   // the coefficient matrix that expresses every beta spatial orbital (as a column) in its underlying scalar basis


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  A basic constructor that sets all the properties.
     * 
     *  @param N_alpha                                  the number of alpha electrons
     *  @param N_beta                                   the number of beta electrons
     *  @param orbital_energies_alpha                   the orbital energies for the alpha-spin-orbitals, sorted by ascending energy
     *  @param orbital_energies_beta                    the orbital energies for the beta-spin-orbitals, sorted by ascending energy
     *  @param C_alpha                                  the coefficient matrix that expresses every alpha spatial orbital (as a column) in its underlying scalar basis
     *  @param C_beta                                   the coefficient matrix that expresses every beta spatial orbital (as a column) in its underlying scalar basis
     */
    UHF(const size_t N_alpha, const size_t N_beta, const VectorX<double>& orbital_energies_alpha, const VectorX<double>& orbital_energies_beta, const TransformationMatrix<Scalar>& C_alpha, const TransformationMatrix<Scalar>& C_beta) :
        N_alpha {N_alpha},
        N_beta {N_beta},
        orbital_energies_alpha {orbital_energies_alpha},
        orbital_energies_beta {orbital_energies_beta},
        C_alpha {C_alpha},
        C_beta {C_beta} {

        // Check for valid arguments.
        const auto K_alpha = C_alpha.dimension();  // number of alpha spatial orbitals
        const auto K_beta = C_beta.dimension();    // number of beta spatial orbitals

        if (N_alpha > K_alpha) {
            throw std::invalid_argument("UHF(const size_t, const size_t, const VectorX<double>&, const VectorX<double>&, const TransformationMatrix<Scalar>&, const TransformationMatrix<Scalar>&): The number of given alpha electrons cannot be larger than the number of alpha spatial orbitals.");
        }

        if (N_beta > K_beta) {
            throw std::invalid_argument("UHF(const size_t, const size_t, const VectorX<double>&, const VectorX<double>&, const TransformationMatrix<Scalar>&, const TransformationMatrix<Scalar>&): The number of given beta electrons cannot be larger than the number of beta spatial orbitals.");
        }

        if (K_alpha != orbital_energies_alpha.size()) {
            throw std::invalid_argument("UHF(const size_t, const size_t, const VectorX<double>&, const VectorX<double>&, const TransformationMatrix<Scalar>&, const TransformationMatrix<Scalar>&): The number of given alpha-spin-orbital energies does not match the number of alpha spin-orbitals.");
        }

        if (K_beta != orbital_energies_beta.size()) {
            throw std::invalid_argument("UHF(const size_t, const size_t, const VectorX<double>&, const VectorX<double>&, const TransformationMatrix<Scalar>&, const TransformationMatrix<Scalar>&): The number of given beta-spin-orbital energies does not match the beta of beta spin-orbitals.");
        }
    }


    /**
     *  A constructor that initializes the orbital energies to zeros.
     * 
     *  @param N_alpha                                  the number of alpha electrons
     *  @param N_beta                                   the number of beta electrons
     *  @param orbital_energies_alpha                   the orbital energies for the alpha-spin-orbitals, sorted by ascending energy
     *  @param orbital_energies_beta                    the orbital energies for the beta-spin-orbitals, sorted by ascending energy
     *  @param C_alpha                                  the coefficient matrix that expresses every alpha spatial orbital (as a column) in its underlying scalar basis
     *  @param C_beta                                   the coefficient matrix that expresses every beta spatial orbital (as a column) in its underlying scalar basis
     */
    UHF(const size_t N_alpha, const size_t N_beta, const TransformationMatrix<Scalar>& C_alpha, const TransformationMatrix<Scalar>& C_beta) :
        UHF(N_alpha, N_beta,
            GQCP::VectorX<double>::Zero(C_alpha.dimension()),
            GQCP::VectorX<double>::Zero(C_beta.dimension()),
            C_alpha, C_beta) {
    }


    /**
     *  Convert an RHF wave function model to an UHF wave function model.
     * 
     *  @param rhf_model            an RHF wave function model
     */
    UHF(const GQCP::QCModel::RHF<Scalar>& rhf_model) :
        UHF(rhf_model.numberOfElectrons(Spin::alpha), rhf_model.numberOfElectrons(Spin::beta),
            rhf_model.orbitalEnergies(), rhf_model.orbitalEnergies(),
            rhf_model.coefficientMatrix(), rhf_model.coefficientMatrix()) {}
};


}  // namespace QCModel
}  // namespace GQCP
