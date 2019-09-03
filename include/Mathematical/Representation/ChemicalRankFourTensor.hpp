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


#include "Basis/TransformationMatrix.hpp"
#include "Mathematical/Representation/SquareRankFourTensor.hpp"

#include <iostream>


namespace GQCP {


/**
 *  An extension of the square rank-4 tensor with methods of quantum chemical context
 *
 *  @tparam _Scalar      the scalar type
 */
template<typename _Scalar>
class ChemicalRankFourTensor: public SquareRankFourTensor<_Scalar> {
public:
    using Scalar = _Scalar;

    using Base = SquareRankFourTensor<Scalar>;
    using Self = ChemicalRankFourTensor<Scalar>;


public:

    /*
     *  CONSTRUCTORS
     */

    using SquareRankFourTensor<Scalar>::SquareRankFourTensor;  // use base constructors


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return the dimension of this matrix representation of the parameters, i.e. the number of orbitals/sites
     */
    size_t dimension() const {
        return static_cast<size_t>(this->Base::dimension(0));  // returns a long
    }

    size_t get_K() const { return this->dimension(0); };



    /**
     *  In-place transform this chemical rank-4 tensor according to a given basis transformation
     *
     *  @param T    the transformation matrix between the old and the new orbital basis, it is used as
     *      b' = b T ,
     *   in which the basis functions are collected as elements of a row vector b
     */
    void basisTransformInPlace(const TransformationMatrix<Scalar>& T) {

        // Since we're only getting T as a matrix, we should make the appropriate tensor to perform contractions
        // For the const argument, we need the const in the template
        //      For more info, see: https://stackoverflow.com/questions/45283468/eigen-const-tensormap
        Eigen::TensorMap<Eigen::Tensor<const Scalar, 2>> T_tensor (T.data(), T.rows(), T.cols());


        // We will have to do four single contractions, so we specify the contraction indices
        // Eigen3 does not document its tensor contraction clearly, so see the accepted answer on stackoverflow (https://stackoverflow.com/a/47558349/7930415):
        //      Eigen3 does not accept a way to specify the output axes: instead, it retains the order from left to right of the axes that survive the contraction.
        //      This means that, in order to get the right ordering of the axes, we will have to swap axes

        // g(T U V W)  T^*(V R) -> a(T U R W) but we get a(T U W R)
        Eigen::array<Eigen::IndexPair<int>, 1> contraction_pair1 = {Eigen::IndexPair<int>(2, 0)};
        Eigen::array<int, 4> shuffle_1 {0, 1, 3, 2};

        // a(T U R W)  T(W S) -> b(T U R S) and we get b(T U R S), so no shuffle is needed
        Eigen::array<Eigen::IndexPair<int>, 1> contraction_pair2 = {Eigen::IndexPair<int>(3, 0)};

        // T(U Q)  b(T U R S) -> c(T Q R S) but we get c(Q T R S)
        Eigen::array<Eigen::IndexPair<int>, 1> contraction_pair3 = {Eigen::IndexPair<int>(0, 1)};
        Eigen::array<int, 4> shuffle_3 {1, 0, 2, 3};

        // T^*(T P)  c(T Q R S) -> g'(P Q R S) and we get g_SO(P Q R S), so no shuffle is needed
        Eigen::array<Eigen::IndexPair<int>, 1> contraction_pair4 = {Eigen::IndexPair<int>(0, 0)};


        // Calculate the contractions. We write this as one large contraction to
        //  1) avoid storing intermediate contractions
        //  2) let Eigen figure out some optimizations
        Self g_transformed = T_tensor.conjugate().contract(T_tensor.contract(this->contract(T_tensor.conjugate(), contraction_pair1).shuffle(shuffle_1).contract(T_tensor, contraction_pair2), contraction_pair3).shuffle(shuffle_3), contraction_pair4);
        (*this) = g_transformed;
    }


    /**
     *  In-place rotate this chemical rank-4 tensor using a unitary transformation matrix
     * 
     *  @param jacobi_rotation_parameters       the Jacobi rotation parameters (p, q, angle) that are used to specify a Jacobi rotation: we use the (cos, sin, -sin, cos) definition for the Jacobi rotation matrix. See transform() for how the transformation matrix between the two bases should be represented
     */
    void basisRotateInPlace(const TransformationMatrix<double>& U) {

        // Check if the given matrix is actually unitary
        if (!U.isUnitary(1.0e-12)) {
            throw std::invalid_argument("ChemicalRankFourTensor::basisRotateInPlace(const TransformationMatrix<Scalar>&): The given transformation matrix is not unitary.");
        }

        this->basisTransformInPlace(U);
    }


    /**
     *  In-place rotate this chemical rank-4 tensor using Jacobi rotation parameters
     * 
     *  @param jacobi_rotation_parameters       the Jacobi rotation parameters (p, q, angle) that are used to specify a Jacobi rotation: we use the (cos, sin, -sin, cos) definition for the Jacobi rotation matrix. See transform() for how the transformation matrix between the two bases should be represented
     */
    void basisRotateInPlace(const JacobiRotationParameters& jacobi_rotation_parameters) {

        /**
         *  While waiting for an analogous Eigen::Tensor Jacobi module, we implement this rotation by constructing a Jacobi rotation matrix and then simply doing a rotation with it
         */

        const auto dim = this->dimension();
        const auto J = TransformationMatrix<double>::FromJacobi(jacobi_rotation_parameters, dim);
        this->basisRotateInPlace(J);
    }
};


}  // namespace GQCP
