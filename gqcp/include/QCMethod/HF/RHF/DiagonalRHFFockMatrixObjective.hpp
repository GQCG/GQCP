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


#include "QCModel/HF/RHF.hpp"


namespace GQCP {


/**
 *  An objective that checks if the RHF Fock matrix is diagonal, i.e. if the RHF parameters represent the canonical RHF coefficients.
 * 
 *  @tparam _Scalar                 the type of scalar used for the expansion coefficients
 */
template <typename _Scalar>
class DiagonalRHFFockMatrixObjective {

public:
    using Scalar = _Scalar;


private:
    double precision;                       // the precision with which the diagonality of the Fock matrix should be checked
    RSQHamiltonian<Scalar> sq_hamiltonian;  // the Hamiltonian expressed in a scalar basis


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param sq_hamiltonian       the Hamiltonian expressed in a scalar basis
     */
    DiagonalRHFFockMatrixObjective(const RSQHamiltonian<Scalar>& sq_hamiltonian, const double precision = 1.0e-08) :
        precision {precision},
        sq_hamiltonian {sq_hamiltonian} {}


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @param rhf_parameters           the RHF model parameters that are considered to be optimized
     * 
     *  @return if this objective is satisfied with the given optimized model parameters
     */
    bool isSatisfiedWith(const QCModel::RHF<Scalar>& rhf_parameters) const {

        // Prepare the calculation of the Fock matrix in the spinor basis
        const auto N_P = rhf_parameters.numberOfElectronPairs();
        const auto C = rhf_parameters.expansion();
        const auto D = QCModel::RHF<Scalar>::calculateScalarBasis1DM(C, 2 * N_P);

        // Calculate the Fock matrix in the orthonormal spinor and check if it is diagonal
        auto F_orthonormal = QCModel::RHF<Scalar>::calculateScalarBasisFockMatrix(D, this->sq_hamiltonian);  // the converged Fock matrix (expressed in the scalar orbital basis)
        F_orthonormal.transform(C);                                                                          // now in the orthonormal spinor basis
        return F_orthonormal.parameters().isDiagonal(this->precision);
    }
};


}  // namespace GQCP
