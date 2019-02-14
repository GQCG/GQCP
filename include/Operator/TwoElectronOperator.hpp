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
#include "math/SquareFourIndexTensor.hpp"
#include "Operator/Operator.hpp"
#include "utilities/miscellaneous.hpp"


namespace GQCP {


/**
 *  A class that represents a two-electron operator in an orbital basis
 *
 *  @tparam Scalar      the scalar type
 */
template<typename Scalar>
class TwoElectronOperator : public SquareFourIndexTensor<Scalar>, public Operator<Scalar> {
public:

    /*
     * CONSTRUCTORS
     */

    /**
     *  Default constructor
     */
    TwoElectronOperator() :
        SquareFourIndexTensor<Scalar>()
    {}


    /**
     *  @param tensor   the explicit matrix representation of the two-electron operator
     */
    explicit TwoElectronOperator(const SquareFourIndexTensor<Scalar>& tensor) :
        SquareFourIndexTensor<Scalar>(tensor)
    {}


    /**
     *  In-place transform the matrix representation of the two-electron operator
     *
     *  @param T    the transformation matrix between the old and the new orbital basis, it is used as
     *      b' = b T ,
     *   in which the basis functions are collected as elements of a row vector b
     */
    void transform(const SquareMatrix<Scalar>& T) override {

        // Since we're only getting T as a matrix, we should make the appropriate tensor to perform contractions
        // For the const Eigen::MatrixXd& argument, we need the const double in the template
        //      For more info, see: https://stackoverflow.com/questions/45283468/eigen-const-tensormap
        Eigen::TensorMap<Eigen::Tensor<const double, 2>> T_tensor (T.data(), T.rows(), T.cols());


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
        //  2) let Eigen3 figure out some optimizations
        Eigen::Tensor<double, 4> g_transformed = T_tensor.conjugate().contract(T_tensor.contract(this->contract(T_tensor.conjugate(), contraction_pair1).shuffle(shuffle_1).contract(T_tensor, contraction_pair2), contraction_pair3).shuffle(shuffle_3), contraction_pair4);

        *this = TwoElectronOperator<Scalar>(g_transformed);
    }


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
        auto J = SquareMatrix<double>(jacobiRotationMatrix(jacobi_rotation_parameters, dim));  // this is sure to return a unitary matrix

        this->Operator<Scalar>::rotate(J);
    }
};



}  // namespace GQCP


#endif  // GQCP_TWOELECTRONOPERATOR_HPP
