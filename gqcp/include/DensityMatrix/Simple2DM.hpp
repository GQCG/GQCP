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


#include "Basis/Transformations/BasisTransformable.hpp"
#include "Basis/Transformations/JacobiRotatable.hpp"
#include "DensityMatrix/DensityMatrixTraits.hpp"
#include "Mathematical/Representation/SquareRankFourTensor.hpp"


namespace GQCP {


/*
 *  MARK: `Simple2DM` implementation
 */

/**
 *  A two-electron density matrix that is described by a single tensor.
 * 
 *  This class is used as a base class for `Orbital2DM` and `G2DM`, since they are both expressed using a single tensor, as opposed to `SpinResolved2DM`, which uses separate alpha- and beta- tensor. The word 'simple' is used here as an antonym for 'compound'.
 * 
 *  @tparam _Scalar                 The scalar type used for a density matrix element: real or complex.
 *  @tparam _DerivedDM              The type of the density matrix that derives from this class, enabling CRTP and compile-time polymorphism.
 */
template <typename _Scalar, typename _DerivedDM>
class Simple2DM:
    public SquareRankFourTensor<_Scalar>,
    public BasisTransformable<_DerivedDM>,
    public JacobiRotatable<_DerivedDM> {
public:
    // The scalar type used for a density matrix element: real or complex.
    using Scalar = _Scalar;

    // The type of the density matrix that derives from this class, enabling CRTP and compile-time polymorphism.
    using DerivedDM = _DerivedDM;

    // The type of 'this'.
    using Self = Simple2DM<Scalar, DerivedDM>;

    // The type of the one-electron density matrix that is naturally related to the derived 2-DM.
    using OneDM = typename DensityMatrixTraits<DerivedDM>::OneDM;

    // The type of transformation that is naturally associated to the derived 2-DM.
    using Transformation = typename DensityMatrixTraits<DerivedDM>::Transformation;


public:
    /*
     *  MARK: Constructors
     */

    // Inherit `SquareRankFourTensor`'s constructors.
    using SquareRankFourTensor<Scalar>::SquareRankFourTensor;


    /*
     *  MARK: General information
     */

    /**
     *  @return The number of orbitals that are related to this 2-DM.
     */
    size_t numberOfOrbitals() const { return this->dimension(); }


    /*
     *  MARK: Contractions
     */

    /**
     *  @return A partial contraction D(p,q) of the 2-DM, where D(p,q) = d(p,q,r,r).
     */
    OneDM reduce() const {
        // TODO: when Eigen3 releases tensor.trace(), use it to implement the reduction

        const auto K = this->numberOfOrbitals();

        OneDM D = OneDM::Zero(K);
        for (size_t p = 0; p < K; p++) {
            for (size_t q = 0; q < K; q++) {

                for (size_t r = 0; r < K; r++) {
                    D(p, q) += this->operator()(p, q, r, r);
                }
            }
        }

        return D;
    }


    /**
     *  @return The trace of the 2-DM, i.e. d(p,p,q,q).
     */
    Scalar trace() const {
        // TODO: when Eigen3 releases tensor.trace(), use it to implement the trace

        const auto K = this->numberOfOrbitals();

        Scalar trace {};
        for (size_t p = 0; p < K; p++) {
            for (size_t q = 0; q < K; q++) {
                trace += this->operator()(p, p, q, q);
            }
        }

        return trace;
    }


    /*
     *  MARK: Conforming to `BasisTransformable`
     */

    /**
     *  Apply the basis transformation and return the resulting two-electron integrals.
     * 
     *  @param T            The basis transformation
     * 
     *  @return The basis-transformed one-electron integrals.
     */
    DerivedDM transformed(const Transformation& T) const override {

        // The transformation formulas for two-electron operators and 2-DMs are similar, but not quite the same. Instead of using T, the transformation formula for the 2-DM uses T_inverse_transpose. See also (https://gqcg-res.github.io/knowdes/spinor-transformations.html).

        // Since we're only getting T as a matrix, we should convert it to an appropriate tensor to perform contractions.
        // Although not a necessity for the einsum implementation, it makes it a lot easier to follow the formulas.
        const GQCP::SquareMatrix<Scalar> T_related = T.matrix().transpose().inverse();
        const GQCP::Tensor<Scalar, 2> T_related_tensor = Eigen::TensorMap<Eigen::Tensor<const Scalar, 2>>(T_related.data(), T_related.rows(), T_related.cols());

        // We calculate the conjugate as a tensor as well.
        const GQCP::Tensor<Scalar, 2> T_related_conjugate = T_related_tensor.conjugate();

        // We will have to do four single contractions
        // d(T U V W)  T^*(V R) -> a(T U R W)
        // a(T U R W)  T(W S) -> b(T U R S)
        const auto temp_1 = this->template einsum<1>("TUVW,VR->TURW", T_related_conjugate).template einsum<1>("TURW,WS->TURS", T_related_tensor);

        // T(U Q)  b(T U R S) -> c(T Q R S)
        const auto temp_2 = T_related_tensor.template einsum<1>("UQ,TURS->TQRS", temp_1);

        // T^*(T P)  c(T Q R S) -> d'(P Q R S)
        const auto transformed = T_related_conjugate.template einsum<1>("TP,TQRS->PQRS", temp_2);

        return DerivedDM {transformed};
    }


    // Allow the `rotate` method from `BasisTransformable`, since there's also a `rotate` from `JacobiRotatable`.
    using BasisTransformable<DerivedDM>::rotate;

    // Allow the `rotated` method from `BasisTransformable`, since there's also a `rotated` from `JacobiRotatable`.
    using BasisTransformable<DerivedDM>::rotated;


    /*
     *  MARK: Conforming to `JacobiRotatable`
     */

    /**
     *  Apply the Jacobi rotation and return the result.
     * 
     *  @param jacobi_rotation          The Jacobi rotation.
     * 
     *  @return The Jacobi-transformed object.
     */
    DerivedDM rotated(const JacobiRotation& jacobi_rotation) const override {

        // While waiting for an analogous Eigen::Tensor Jacobi module, we implement this rotation by constructing a Jacobi rotation matrix and then simply doing a rotation with it.
        const auto J = Transformation::FromJacobi(jacobi_rotation, this->numberOfOrbitals());
        return this->rotated(J);
    }

    // Allow the `rotate` method from `JacobiRotatable`, since there's also a `rotate` from `BasisTransformable`.
    using JacobiRotatable<DerivedDM>::rotate;
};


}  // namespace GQCP
