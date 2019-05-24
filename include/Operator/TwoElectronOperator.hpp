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
#ifndef GQCP_TWOELECTRONOPERATOR_HPP
#define GQCP_TWOELECTRONOPERATOR_HPP


#include "JacobiRotationParameters.hpp"
#include "math/ChemicalRankFourTensor.hpp"
#include "Operator/OneElectronOperator.hpp"
#include "utilities/miscellaneous.hpp"


namespace GQCP {


/**
 *  A class that represents a two-electron operator in an orbital basis
 *
 *  @tparam _Scalar     the scalar type
 */
template<typename _Scalar>
class TwoElectronOperator : public ChemicalRankFourTensor<_Scalar>, public Operator<TwoElectronOperator<_Scalar>> {
public:

    using Scalar = _Scalar;

    using BaseRepresentation = SquareMatrix<Scalar>;
    using Self = TwoElectronOperator<Scalar>;


public:

    /*
     *  CONSTRUCTORS
     */

    using ChemicalRankFourTensor<Scalar>::ChemicalRankFourTensor;  // use base constructors



    /*
     *  PUBLIC METHODS
     */


    using Operator<TwoElectronOperator<Scalar>>::rotate;  // bring over rotate() from the base class


    /**
     *  In-place rotate the matrix representation of the two-electron operator using a unitary Jacobi rotation matrix constructed from the Jacobi rotation parameters. Note that this function is only available for real (double) matrix representations
     *
     *  @param jacobi_rotation_parameters       the Jacobi rotation parameters (p, q, angle) that are used to specify a Jacobi rotation: we use the (cos, sin, -sin, cos) definition for the Jacobi rotation matrix. See transform() for how the transformation matrix between the two bases should be represented
     */
    template<typename Z = Scalar>
    enable_if_t<std::is_same<Z, double>::value> rotate(const JacobiRotationParameters& jacobi_rotation_parameters) {

        /**
         *  While waiting for an analogous Eigen::Tensor Jacobi module, we implement this rotation by constructing a Jacobi rotation matrix and then doing a rotation with it
         */

        auto dim = static_cast<size_t>(this->dimension(0));  // .dimension() returns a long
        auto J = SquareMatrix<double>::FromJacobi(jacobi_rotation_parameters, dim);  // this is sure to return a unitary matrix

        this->rotate(J);
    }


    /**
     *  @return the two-electron integrals that can be evaluated through a one electron mode as a one-electron operator
     */
    OneElectronOperator<Scalar> effectiveOneElectronPartition() const {

        auto K = this->dimension(0);

        OneElectronOperator<Scalar> k = OneElectronOperator<Scalar>::Zero(K, K);

        for (size_t p = 0; p < K; p++) {
            for (size_t q = 0; q < K; q++) {
                for (size_t r = 0; r < K; r++) {
                    k(p, q) -= 0.5 * this->operator()(p, r, r, q);
                }
            }
        }

        return k;
    }
};



}  // namespace GQCP


#endif  // GQCP_TWOELECTRONOPERATOR_HPP
