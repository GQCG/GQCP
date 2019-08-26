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

#include "Mathematical/ChemicalRankFourTensor.hpp"
#include "Operator/SecondQuantized/SQOneElectronOperator.hpp"
#include "OrbitalOptimization/JacobiRotationParameters.hpp"
#include "Utilities/miscellaneous.hpp"

#include <array>


namespace GQCP {


/**
 *  A class that represents a second-quantized two-electron operator: it holds the matrix representation of its parameters, which are (usually) integrals over first-quantized operators
 *
 *  @tparam _Scalar             the scalar type, i.e. the scalar representation of one of the parameters
 *  @tparam _Components         the number of components of the second-quantized operator
 */
template <typename _Scalar, size_t _Components>
class SQTwoElectronOperator {
public:

    using Scalar = _Scalar;
    static constexpr auto Components = _Components;


private:
    std::array<ChemicalRankFourTensor<Scalar>, Components> G;  // all the matrix representations of the parameters (integrals) of the different components of this second-quantized operator

public:

    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param G            all the matrix representations of the parameters (integrals) of the different components of the second-quantized operator
     */
    SQTwoElectronOperator(const std::array<ChemicalRankFourTensor<Scalar>, Components>& G) : 
        G (G)
    {
        // Check if the given matrix representations have the same dimensions
        const auto dimension = this->G[0].dimension();

        for (size_t i = 1; i < Components; i++) {
            if (dimension != this->G[i].dimension()) {
                throw std::invalid_argument("SQOneElectronOperator(const std::array<ChemicalMatrix<Scalar>, Components>&): The given matrix representations did not have the same dimensions.");
            }
        }
    }


    /**
     *  Construct a two-electron operator with zero parameters
     * 
     *  @param dim          the dimension of the matrix representation of the parameters, i.e. the number of orbitals/sites
     */
    SQTwoElectronOperator(const size_t dim) {
        for (size_t i = 0; i < Components; i++) {
            this->G[i] = ChemicalRankFourTensor<Scalar>(dim);
        }
    }



    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return the dimension of the matrix representation of the parameters, i.e. the number of orbitals/sites
     */
    size_t dimension() const {
        return this->G[0].dimension();
    }

    size_t get_dim() const {
        return this->dimension();
    }

    size_t get_K() const {
        return this->dimension();
    }


    /**
     *  @return all the matrix representations of the parameters (integrals) of the different components of this second-quantized operator
     */
    const std::array<ChemicalRankFourTensor<Scalar>, Components> allParameters() const {
        return this->G;
    }


    /**
     *  @param i            the index of the component
     * 
     *  @return a read-only matrix representation of the parameters (integrals) of one of the the different components of this second-quantized operator
     */
    const ChemicalRankFourTensor<Scalar>& parameters(const size_t i = 0) const {
        return this->G[i];
    }


    /**
     *  @param i            the index of the component
     * 
     *  @return a writable the matrix representation of the parameters (integrals) of one of the the different components of this second-quantized operator
     */
    ChemicalRankFourTensor<Scalar>& parameters(const size_t i = 0) {
        return this->G[i];
    }


    /**
     *  @tparam TransformationScalar        the scalar type of the transformation matrix
     * 
     *  @param T                            the transformation matrix
     * 
     *  @return a second-quantized operator whose matrix representations have been transformed according to the given transformation matrix
     */
    template <typename TransformationScalar = Scalar>
    auto transform(const SquareMatrix<TransformationScalar>& T) const -> SQTwoElectronOperator<product_t<Scalar, TransformationScalar>, Components> const {

        using ResultScalar = product_t<Scalar, TransformationScalar>;

        // Transform the matrix representations of the components
        auto G_copy = this->allParameters();
        for (size_t i = 0; i < Components; i++) {
            G_copy[i] = this->parameters(i).basisTransform(T);
        }

        return SQTwoElectronOperator<ResultScalar, Components>(G_copy);
    }


    /**
     *  @tparam TransformationScalar        the scalar type of the transformation matrix
     * 
     *  @param U                            the (unitary) rotation matrix
     * 
     *  @return a second-quantized operator whose matrix representations have been rotated according to the given rotation matrix
     */
    template <typename TransformationScalar = Scalar>
    auto rotate(const SquareMatrix<TransformationScalar>& U) const -> SQTwoElectronOperator<product_t<Scalar, TransformationScalar>, Components> const {

        // Check if the given matrix is actually unitary
        if (!U.isUnitary(1.0e-12)) {
            throw std::invalid_argument("SQTwoElectronOperator::rotate(const SquareMatrix<TransformationScalar>&): The given transformation matrix is not unitary.");
        }

        return this->transform(U);
    }


    /**
     *  In-place rotate the matrix representation of the two-electron operator using a unitary Jacobi rotation matrix constructed from the Jacobi rotation parameters. Note that this function is only available for real (double) matrix representations
     *
     *  @param jacobi_rotation_parameters       the Jacobi rotation parameters (p, q, angle) that are used to specify a Jacobi rotation: we use the (cos, sin, -sin, cos) definition for the Jacobi rotation matrix. See transform() for how the transformation matrix between the two bases should be represented
     */
    template<typename Z = Scalar>
    enable_if_t<std::is_same<Z, double>::value,
    SQTwoElectronOperator<product_t<Z, double>, Components>> rotate(const JacobiRotationParameters& jacobi_rotation_parameters) const {

        using ResultScalar = product_t<Scalar, double>;

        // Transform the matrix representations of the components
        auto G_copy = this->allParameters();
        for (size_t i = 0; i < Components; i++) {
            G_copy[i] = this->parameters(i).basisRotate(jacobi_rotation_parameters);
        }

        return SQTwoElectronOperator<ResultScalar, Components>(G_copy);
    }


    /**
     *  @return the one-electron operator that is the difference between a two-electron operator (e_pqrs) and a product of one-electron operators (E_pq E_rs)
     */
    SQOneElectronOperator<Scalar, Components> effectiveOneElectronPartition() const {

        // Initialize a zero operator
        const auto K = this->dimension();  // number of orbitals
        SQOneElectronOperator<Scalar, Components> F (K);


        // Use a formula to set the parameters
        for (size_t i = 0; i < Components; i++) {
            for (size_t p = 0; p < K; p++) {
                for (size_t q = 0; q < K; q++) {
                    for (size_t r = 0; r < K; r++) {
                        F.parameters(i)(p,q) -= 0.5 * this->parameters(i)(p, r, r, q);
                    }
                }
            }
        }

        return F;
    }
};



/*
 *  CONVENIENCE ALIASES
 */
template <typename Scalar>
using ScalarSQTwoElectronOperator = SQTwoElectronOperator<Scalar, 1>;

template <typename Scalar>
using VectorSQTwoElectronOperator = SQTwoElectronOperator<Scalar, 3>;


/*
 *  OPERATORS
 */

/**
 *  Add two two-electron operators by adding their parameters
 * 
 *  @tparam LHSScalar           the scalar type of the left-hand side
 *  @tparam RHSScalar           the scalar type of the right-hand side
 *  @tparam Components          the number of components of the two-electron operators
 * 
 *  @param lhs                  the left-hand side
 *  @param rhs                  the right-hand side
 */
template <typename LHSScalar, typename RHSScalar, size_t Components>
auto operator+(const SQTwoElectronOperator<LHSScalar, Components>& lhs, const SQTwoElectronOperator<RHSScalar, Components>& rhs) -> SQTwoElectronOperator<sum_t<LHSScalar, RHSScalar>, Components> {

    using ResultScalar = sum_t<LHSScalar, RHSScalar>;

    auto G_sum = lhs.allParameters();
    for (size_t i = 0; i < Components; i++) {
        G_sum[i] += rhs.parameters(i);
    }

    return SQTwoElectronOperator<ResultScalar, Components>(G_sum);
}


/**
 *  Multiply a two-electron operator with a scalar
 * 
 *  @tparam Scalar              the scalar type of the scalar
 *  @tparam OperatorScalar      the scalar type of the operator
 * 
 *  @tparam scalar              the scalar of the scalar multiplication
 *  @tparam op                  the two-electron operator
 */
template <typename Scalar, typename OperatorScalar, size_t Components>
auto operator*(const Scalar& scalar, const SQTwoElectronOperator<OperatorScalar, Components>& op) -> SQTwoElectronOperator<product_t<Scalar, OperatorScalar>, Components> {

    using ResultScalar = product_t<Scalar, OperatorScalar>;

    auto G = op.allParameters();
    for (size_t i = 0; i < Components; i++) {
        G[i] *= scalar;
    }

    return SQTwoElectronOperator<ResultScalar, Components>(G);
}


}  // namespace GQCP
