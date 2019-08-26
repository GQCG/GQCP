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


#include "Mathematical/ChemicalMatrix.hpp"
#include "Mathematical/ScalarFunction.hpp"
#include "OrbitalOptimization/JacobiRotationParameters.hpp"
#include "typedefs.hpp"

#include <array>


namespace GQCP {


/**
 *  A class that represents a second-quantized one-electron operator: it holds the matrix representation of its parameters, which are (usually) integrals over first-quantized operators
 *
 *  @tparam _Scalar             the scalar type, i.e. the scalar representation of one of the parameters
 *  @tparam _Components         the number of components of the second-quantized operator
 */
template <typename _Scalar, size_t _Components>
// class SQOneElectronOperator : public BaseSQOperator<SQOneElectronOperator<_Scalar>> {
class SQOneElectronOperator {
public:

    using Scalar = _Scalar;
    static constexpr auto Components = _Components;


private:
    std::array<ChemicalMatrix<Scalar>, Components> F;  // all the matrix representations of the parameters (integrals) of the different components of this second-quantized operator


public:

    /*
     *  CONSTRUCTORS
     */

    /**
     *  
     * 
     *  @param F            all the matrix representations of the parameters (integrals) of the different components of the second-quantized operator
     */
    SQOneElectronOperator(const std::array<ChemicalMatrix<Scalar>, Components>& F) : 
        F (F)
    {
        // Check if the given matrix representations have the same dimensions
        const auto dimension = this->F[0].dimension();

        for (size_t i = 1; i < Components; i++) {
            if (dimension != this->F[i].dimension()) {
                throw std::invalid_argument("SQOneElectronOperator(const std::array<ChemicalMatrix<Scalar>, Components>&): The given matrix representations did not have the same dimensions.");
            }
        }
    }


    /**
     *  Construct a one-electron operator with parameters that are zero
     * 
     *  @param dim          the dimension of the matrix representation of the parameters, i.e. the number of orbitals/sites
     */
    SQOneElectronOperator(const size_t dim) {
        for (size_t i = 0; i < Components; i++) {
            this->F[i] = ChemicalMatrix<Scalar>(dim);
        }
    }


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return the dimension of the matrix representation of the parameters, i.e. the number of orbitals/sites
     */
    size_t dimension() const {
        return this->F[0].dimension();
    }

    size_t get_dim() const {
        return this->dimension();
    }


    /**
     *  @return all the matrix representations of the parameters (integrals) of the different components of this second-quantized operator
     */
    const std::array<ChemicalMatrix<Scalar>, Components> allParameters() const {
        return this->F;
    }


    /**
     *  @param i            the index of the component
     * 
     *  @return a read-only the matrix representation of the parameters (integrals) of one of the the different components of this second-quantized operator
     */
    const ChemicalMatrix<Scalar>& parameters(const size_t i = 0) const {
        return this->F[i];
    }


    /**
     *  @param i            the index of the component
     * 
     *  @return a writable matrix representation of the parameters (integrals) of one of the the different components of this second-quantized operator
     */
    ChemicalMatrix<Scalar>& parameters(const size_t i = 0) {
        return this->F[i];
    }


    /**
     *  @tparam TransformationScalar        the scalar type of the transformation matrix
     * 
     *  @param T                            the transformation matrix
     * 
     *  @return a second-quantized operator whose matrix representations have been transformed according to the given transformation matrix
     */
    template <typename TransformationScalar = Scalar>
    auto transform(const SquareMatrix<TransformationScalar>& T) const -> SQOneElectronOperator<product_t<Scalar, TransformationScalar>, Components> const {

        using ResultScalar = product_t<Scalar, TransformationScalar>;

        // Transform the matrix representations of the components
        auto F_copy = this->allParameters();
        for (size_t i = 0; i < Components; i++) {
            F_copy[i] = this->parameters(i).basisTransform(T);
        }

        return SQOneElectronOperator<ResultScalar, Components>(F_copy);
    }


    /**
     *  @tparam TransformationScalar        the scalar type of the transformation matrix
     * 
     *  @param U                            the (unitary) rotation matrix
     * 
     *  @return a second-quantized operator whose matrix representations have been rotated according to the given rotation matrix
     */
    template <typename TransformationScalar = Scalar>
    auto rotate(const SquareMatrix<TransformationScalar>& U) const -> SQOneElectronOperator<product_t<Scalar, TransformationScalar>, Components> const {

        // Check if the given matrix is actually unitary
        if (!U.isUnitary(1.0e-12)) {
            throw std::invalid_argument("SQOneElectronOperator::rotate(const SquareMatrix<TransformationScalar>&): The given transformation matrix is not unitary.");
        }

        return this->transform(U);
    }


    /**
     *  Rotate the parameters of this one-electron operator using a unitary Jacobi rotation matrix constructed from the Jacobi rotation parameters. Note that this function is only available for real (double) matrix representations
     *
     *  @param jacobi_rotation_parameters       the Jacobi rotation parameters (p, q, angle) that are used to specify a Jacobi rotation: we use the (cos, sin, -sin, cos) definition for the Jacobi rotation matrix. See transform() for how the transformation matrix between the two bases should be represented
     */
    template<typename Z = Scalar>
    enable_if_t<std::is_same<Z, double>::value,
    SQOneElectronOperator<product_t<Z, double>, Components>> rotate(const JacobiRotationParameters& jacobi_rotation_parameters) const {

        using ResultScalar = product_t<Scalar, double>;

        // Transform the matrix representations of the components
        auto F_copy = this->allParameters();
        for (size_t i = 0; i < Components; i++) {
            F_copy[i] = this->parameters(i).basisRotate(jacobi_rotation_parameters);
        }

        return SQOneElectronOperator<ResultScalar, Components>(F_copy);
    }


    /**
     *  @param x        the vector/point at which the scalar functions should be evaluated
     *
     *  @return a one-electron operator corresponding to the evaluated scalar functions
     *
     *  Note that this function is only available for SQOneElectronOperators whose Scalar is a derived class of ScalarFunction
     */
    template <typename Z = Scalar>
    enable_if_t<std::is_base_of<ScalarFunction<typename Z::Valued, typename Z::Scalar, Z::Cols>, Z>::value,
    SQOneElectronOperator<typename Z::Valued, Components>> evaluate(const Vector<typename Z::Scalar, Z::Cols>& x) const {

        // Initialize the results
        SQOneElectronOperator<typename Z::Value, Components> F_evaluated (this->dimension());

        // Evaluate all components at the given x
        for (size_t i = 0; i < Components; i++) {
            for (size_t m = 0; m < this->rows(); m++) {
                for (size_t n = 0; n < this->cols(); n++) {
                    const auto F_i_mn = this->F[i](m,n);  // (m,n)-th element of the i-th component
                    F_evaluated[i](m,n) = F_i_mn.operator()(x);  // evaluate the ScalarFunction
                }
            }
        }

        return F_evaluated;
    }
};



/*
 *  CONVENIENCE ALIASES
 */
template <typename Scalar>
using ScalarSQOneElectronOperator = SQOneElectronOperator<Scalar, 1>;

template <typename Scalar>
using VectorSQOneElectronOperator = SQOneElectronOperator<Scalar, 3>;


/*
 *  OPERATORS
 */

/**
 *  Add two one-electron operators by adding their parameters
 * 
 *  @tparam LHSScalar           the scalar type of the left-hand side
 *  @tparam RHSScalar           the scalar type of the right-hand side
 *  @tparam Components          the number of components of the one-electron operators
 * 
 *  @param lhs                  the left-hand side
 *  @param rhs                  the right-hand side
 */
template <typename LHSScalar, typename RHSScalar, size_t Components>
auto operator+(const SQOneElectronOperator<LHSScalar, Components>& lhs, const SQOneElectronOperator<RHSScalar, Components>& rhs) -> SQOneElectronOperator<sum_t<LHSScalar, RHSScalar>, Components> {

    using ResultScalar = sum_t<LHSScalar, RHSScalar>;

    auto F_sum = lhs.allParameters();
    for (size_t i = 0; i < Components; i++) {
        F_sum[i] += rhs.parameters(i);
    }

    return SQOneElectronOperator<ResultScalar, Components>(F_sum);
}


/**
 *  Multiply a one-electron operator with a scalar
 * 
 *  @tparam Scalar              the scalar type of the scalar
 *  @tparam OperatorScalar      the scalar type of the operator
 * 
 *  @tparam scalar              the scalar of the scalar multiplication
 *  @tparam op                  the one-electron operator
 */
template <typename Scalar, typename OperatorScalar, size_t Components>
auto operator*(const Scalar& scalar, const SQOneElectronOperator<OperatorScalar, Components>& op) -> SQOneElectronOperator<product_t<Scalar, OperatorScalar>, Components> {

    using ResultScalar = product_t<Scalar, OperatorScalar>;

    auto F = op.allParameters();
    for (size_t i = 0; i < Components; i++) {
        F[i] *= scalar;
    }

    return SQOneElectronOperator<ResultScalar, Components>(F);
}


}  // namespace GQCP
