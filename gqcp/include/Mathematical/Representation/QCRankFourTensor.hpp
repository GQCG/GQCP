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


#include "Basis/Transformations/TransformationMatrix.hpp"
#include "Mathematical/Representation/SquareRankFourTensor.hpp"

#include <iostream>


namespace GQCP {


/**
 *  An extension of the square rank-4 tensor with methods of quantum chemical context for two-electron integrals.
 *
 *  @tparam _Scalar      the scalar type
 */
template <typename _Scalar>
class QCRankFourTensor: public SquareRankFourTensor<_Scalar> {
public:
    using Scalar = _Scalar;

    using Base = SquareRankFourTensor<Scalar>;
    using Self = QCRankFourTensor<Scalar>;


private:
    bool is_antisymmetrized = false;                   // if the two-electron integrals are modified to obey antisymmetry w.r.t. creation and annihilation indices
    bool is_expressed_using_chemists_notation = true;  // if the two-electron integrals are expressed as g_PQRS or (PQ|RS)


public:
    /*
     *  CONSTRUCTORS
     */

    using SquareRankFourTensor<Scalar>::SquareRankFourTensor;  // use base constructors


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return an antisymmetrized version of these two-electron integrals
     * 
     *  @note If these integrals are expressed using
     *      - chemist's notation g_{PQRS}, return g_{PQRS} - g_{PSRQ}
     *      - physicist's notation <PQ|RS>, return <PQ||RS> = <PQ|RS> - <PQ|SR>
     */
    Self antisymmetrized() const {

        // Attempt to modify a copy of these integrals if they haven't been antisymmetrized already.
        auto copy = *this;
        if (!(this->isAntisymmetrized())) {

            if (this->isExpressedUsingChemistsNotation()) {
                Eigen::array<int, 4> shuffle_indices {0, 3, 2, 1};
                copy -= this->shuffle(shuffle_indices);
            }

            else {  // expressed using physicist's notation

                Eigen::array<int, 4> shuffle_indices {0, 1, 3, 2};
                copy -= this->shuffle(shuffle_indices);
            }

            copy.is_antisymmetrized = true;
        }

        return copy;
    }


    /**
     *  In-place antisymmetrize these two-electron integrals.
     * 
     *  @note If these integrals are expressed using
     *      - chemist's notation g_{PQRS}, they are modified to g_{PQRS} - g_{PSRQ}
     *      - physicist's notation <PQ|RS>, they are modified to <PQ||RS> = <PQ|RS> - <PQ|SR>
     */
    void antisymmetrize() { *this = this->antisymmetrized(); }


    /**
     *  In-place rotate this quantum chemical rank-4 tensor using a unitary transformation matrix.
     * 
     *  @param U            the unitary transformation matrix
     */
    void basisRotate(const TransformationMatrix<double>& U) {

        // Check if the given matrix is actually unitary
        if (!U.isUnitary(1.0e-12)) {
            throw std::invalid_argument("QCRankFourTensor::basisRotate(const TransformationMatrix<Scalar>&): The given transformation matrix is not unitary.");
        }

        this->basisTransform(U);
    }


    /**
     *  In-place rotate this chemical rank-4 tensor using Jacobi rotation parameters.
     * 
     *  @param jacobi_rotation_parameters       the Jacobi rotation parameters (p, q, angle) that are used to specify a Jacobi rotation: we use the (cos, sin, -sin, cos) definition for the Jacobi rotation matrix
     */
    void basisRotate(const JacobiRotationParameters& jacobi_rotation_parameters) {

        /**
         *  While waiting for an analogous Eigen::Tensor Jacobi module, we implement this rotation by constructing a Jacobi rotation matrix and then simply doing a rotation with it
         */
        const auto J = TransformationMatrix<double>::FromJacobi(jacobi_rotation_parameters, this->numberOfOrbitals());
        this->basisRotate(J);
    }


    /**
     *  In-place transform this quantum chemical rank-4 tensor according to a given basis transformation.
     *
     *  @param T    the transformation matrix between the old and the new orbital basis
     */
    void basisTransform(const TransformationMatrix<Scalar>& T) {

        // Since we're only getting T as a matrix, we should convert it to an appropriate tensor to perform contractions.
        Eigen::TensorMap<Eigen::Tensor<const Scalar, 2>> T_tensor {T.data(), T.rows(), T.cols()};  // since T is const, we need const in the template (https://stackoverflow.com/questions/45283468/eigen-const-tensormap)


        // We will have to do four single contractions, so we'll have to specify the contraction indices.
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


        // Calculate the contractions. We write this as one chain of contractions to
        //      1) avoid storing intermediate contractions;
        //      2) let Eigen figure out some optimizations.
        Self g_transformed = T_tensor.conjugate().contract(
            T_tensor.contract(
                        this->contract(T_tensor.conjugate(), contraction_pair1).shuffle(shuffle_1)  // the 'inner' contraction, the first one
                            .contract(T_tensor, contraction_pair2),
                        contraction_pair3)
                .shuffle(shuffle_3),
            contraction_pair4);
        (*this) = g_transformed;
    }


    /**
     *  @return the integrals changed to chemist's notation (from physicist's notation).
     */
    Self convertedToChemistsNotation() const {

        // Attempt to modify a copy if these integrals are expressed in physicist's notation.
        auto copy = *this;
        if (this->isExpressedUsingPhysicistsNotation()) {

            Eigen::array<int, 4> shuffle_indices {0, 2, 1, 3};
            copy = QCRankFourTensor<double>(this->shuffle(shuffle_indices));

            copy.is_expressed_using_chemists_notation = true;
        }
        return copy;
    }


    /**
     *  @return the integrals changed to physicist's notation (from chemist's notation).
     */
    Self convertedToPhysicistsNotation() const {

        // Attempt to modify a copy if these integrals are expressed in chemist's notation.
        auto copy = *this;
        if (this->isExpressedUsingChemistsNotation()) {

            Eigen::array<int, 4> shuffle_indices {0, 2, 1, 3};
            copy = QCRankFourTensor<double>(this->shuffle(shuffle_indices));

            copy.is_expressed_using_chemists_notation = false;
        }
        return copy;
    }


    /**
     *  In-place change the integrals to chemist's notation (from physicist's notation).
     */
    void convertToChemistsNotation() { *this = this->convertedToChemistsNotation(); }

    /**
     *  In-place change the integrals to physicist's notation (from chemist's notation).
     */
    void convertToPhysicistsNotation() { *this = this->convertedToPhysicistsNotation(); }

    /**
     *  @return if these two-electron integrals are considered to be antisymmetrized.
     * 
     *  @note If so, these integrals represent:
     *      - if they are expressed using chemist's notation:       g_{PQRS} - g_{PSRQ}, i.e. they are antisymmetric upon interchanging the indices PR or QS
     *      - if they are expressed using physicist's notation:     <PQ|RS> - <PQ|SR>, i.e. they are antisymmetric upon interchanging the indices PQ or RS
     */
    bool isAntisymmetrized() const { return this->is_antisymmetrized; }

    /**
     *  @return if these two-electron integrals are expressed using chemist's notation g_{PQRS}, i.e. (PQ|RS)
     */
    bool isExpressedUsingChemistsNotation() const { return this->is_expressed_using_chemists_notation; }

    /**
     *  @return if these two-electron integrals are expressed using physicist's notation <PQ|RS>
     */
    bool isExpressedUsingPhysicistsNotation() const { return !(this->isExpressedUsingChemistsNotation()); }

    /**
     *  @return the number of orbitals (spinors or spin-orbitals, depending on the context) this quantum chemical rank-four tensor is associated to
     */
    size_t numberOfOrbitals() const { return this->dimension(); /* the dimension of the underlying square rank-four tensor */ }
};


}  // namespace GQCP
