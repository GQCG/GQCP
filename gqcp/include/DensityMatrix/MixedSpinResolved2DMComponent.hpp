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


#include "Basis/Transformations/UTransformationComponent.hpp"
#include "Mathematical/Representation/SquareRankFourTensor.hpp"
#include "QuantumChemical/Spin.hpp"


namespace GQCP {


/**
 *  One of the mixed (i.e. alpha-beta or beta-alpha) spin components of a spin-resolved 2-DM.
 * 
 *  @tparam _Scalar                 The scalar type used for a density matrix element: real or complex.
 */
template <typename _Scalar>
class MixedSpinResolved2DMComponent:
    public SquareRankFourTensor<_Scalar> {
public:
    // The scalar type used for a density matrix element: real or complex.
    using Scalar = _Scalar;

    // The type of 'this'.
    using Self = MixedSpinResolved2DMComponent<Scalar>;


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
     *  @return The trace of the 2-DM, i.e. d(p,p,q,q).
     */
    Scalar trace() const {

        // FIXME: This is duplicate code from `Simple2DM`. It would be double work to introduce a temporary intermediate class before resolving issue #559 (https://github.com/GQCG/GQCP/issues/559), which is why we chose to just copy-paste the implementation.

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
     *  MARK: Basis transformations
     */

    /**
     *  Apply the basis transformation for the spin component sigma, and return the resulting 2-DM.
     * 
     *  @param T                            The basis transformation.
     *  @param sigma                        Alpha indicates a transformation of the first two axes, beta indicates a transformation of the second two axes.
     * 
     *  @return The basis-transformed 2-DM.
     * 
     *  @note We apologize for this half-baked API. It is currently present in the code, while issue #559 (https://github.com/GQCG/GQCP/issues/688) is being implemented.
     */
    Self transformed(const UTransformationComponent<Scalar>& T, const Spin sigma) const {

        // The transformation formulas for two-electron operators and 2-DMs are similar, but not quite the same. Instead of using T, the transformation formula for the 2-DM uses T_inverse_transpose. See also (https://gqcg-res.github.io/knowdes/spinor-transformations.html).

        // Since we're only getting T as a matrix, we should convert it to an appropriate tensor to perform contractions.
        // Although not a necessity for the einsum implementation, it makes it a lot easier to follow the formulas.
        const GQCP::SquareMatrix<Scalar> T_related = T.matrix().transpose().inverse();
        const GQCP::Tensor<Scalar, 2> T_related_tensor = Eigen::TensorMap<Eigen::Tensor<const Scalar, 2>>(T_related.data(), T_related.rows(), T_related.cols());

        // We calculate the conjugate as a tensor as well.
        const GQCP::Tensor<Scalar, 2> T_related_conjugate = T_related_tensor.conjugate();


        // Depending on the given spin-component, we should either transform the first two, or the second two axes.
        switch (sigma) {
        case Spin::alpha: {
            const auto temp = T_related_tensor.template einsum<1>("UQ,TUVW->TQVW", *this);
            const auto transformed = T_related_conjugate.template einsum<1>("TP,TQVW->PQVW", temp);
            return Self {transformed};
            break;
        }

        case Spin::beta: {
            const auto temp = this->template einsum<1>("PQVW, WS->PQVS", T_related_tensor);
            const auto transformed = temp.template einsum<1>("PQVS, VR->PQRS", T_related_conjugate);
            return Self {transformed};
            break;
        }
        }
    }


    /**
     *  In-place apply the basis transformation for the spin component sigma.
     * 
     *  @param T                    The basis transformation.
     *  @param sigma                Alpha indicates a transformation of the first two axes, beta indicates a transformation of the second two axes.
     * 
     *  @note We apologize for this half-baked API. It is currently present in the code, while issue #559 (https://github.com/GQCG/GQCP/issues/688) is being implemented.
     */
    void transform(const UTransformationComponent<Scalar>& T, const Spin sigma) {

        auto result = this->transformed(T, sigma);
        *this = result;
    }


    /**
     *  Apply the basis rotation for the spin component sigma, and return the resulting two-electron integrals.
     * 
     *  @param jacobi_rotation              The Jacobi rotation.
     *  @param sigma                        Alpha indicates a transformation of the first two axes, beta indicates a transformation of the second two axes.
     * 
     *  @return The basis-rotated two-electron operator.
     * 
     *  @note We apologize for this half-baked API. It is currently present in the code, while issue #559 (https://github.com/GQCG/GQCP/issues/688) is being implemented.
     */
    Self rotated(const JacobiRotation& jacobi_rotation, const Spin sigma) const {

        const auto J = UTransformationComponent<Scalar>::FromJacobi(jacobi_rotation, this->numberOfOrbitals());
        return this->transformed(J, sigma);
    }


    /**
     *  In-place apply the basis rotation for the spin component sigma, and return the resulting two-electron integrals.
     * 
     *  @param jacobi_rotation              The Jacobi rotation.
     *  @param sigma                        Alpha indicates a transformation of the first two axes, beta indicates a transformation of the second two axes.
     * 
     *  @note We apologize for this half-baked API. It is currently present in the code, while issue #559 (https://github.com/GQCG/GQCP/issues/688) is being implemented.
     */
    void rotate(const JacobiRotation& jacobi_rotation, const Spin sigma) {

        auto result = this->rotated(jacobi_rotation, sigma);
        *this = result;
    }
};


}  // namespace GQCP
