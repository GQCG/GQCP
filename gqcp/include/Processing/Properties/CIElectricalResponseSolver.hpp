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

#include "Processing/Properties/BaseElectricalResponseSolver.hpp"
#include "QCModel/CI/LinearExpansion.hpp"


namespace GQCP {


/**
 *  A class whose instances can solve the response equations for CI.
 */
template <typename _ONVBasis>
class CIElectricalResponseSolver: public BaseElectricalResponseSolver {
public:
    using ONVBasis = _ONVBasis;


private:
    LinearExpansion<ONVBasis> linear_expansion;  // the CI linear expansion


public:
    // CONSTRUCTORS

    /**
     *  @param linear_expansion             the CI linear expansion
     */
    CIElectricalResponseSolver(const LinearExpansion<ONVBasis>& linear_expansion) :
        linear_expansion {linear_expansion} {}


    // PUBLIC OVERRIDDEN METHODS

    /**
     *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal orbital basis
     * 
     *  @return the parameter response constant (k_p), i.e. the second-order parameter partial derivative of the CI energy function
     */
    SquareMatrix<double> calculateParameterResponseConstant(const RSQHamiltonian<double>& sq_hamiltonian) const override {

        // k_p for CI models is just the electronic Hamiltonian evaluated in the ONV basis
        return 2 * this->linear_expansion.onvBasis().evaluateOperatorDense(sq_hamiltonian, true);  // true: need to calculate diagonal values as well
    }

    /**
     *  @param dipole_op                    the dipole integrals in an orthonormal orbital basis
     * 
     *  @return the parameter response force (F_p), i.e. the first-order parameter partial derivative of the perturbation derivative of the CI energy function
     */
    Matrix<double, Dynamic, 3> calculateParameterResponseForce(const VectorRSQOneElectronOperator<double>& dipole_op) const override {

        // Prepare some variables.
        const auto dim = this->linear_expansion.onvBasis().dimension();
        const auto& c = this->linear_expansion.coefficients();

        // F_p for CI is related to the matrix-vector product of the total dipole operator and the CI expansion.
        Matrix<double, Dynamic, 3> F_p = Matrix<double, Dynamic, 3>::Zero(dim, 3);
        for (size_t m = 0; m < 3; m++) {  // m loops over the components of the electrical dipole
            const auto& mu_m = dipole_op[m];

            // Calculate the matrix-vector product of the dipole operator and the CI expansion.
            const auto diagonal = this->linear_expansion.onvBasis().evaluateOperatorDiagonal(mu_m);
            F_p.col(m) = -2 * this->linear_expansion.onvBasis().evaluateOperatorMatrixVectorProduct(mu_m, c, diagonal);
        }

        return F_p;
    }
};


}  // namespace GQCP
