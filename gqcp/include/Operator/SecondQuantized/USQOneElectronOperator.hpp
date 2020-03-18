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

#include "Basis/SpinorBasis/JacobiRotationParameters.hpp"
#include "Basis/SpinorBasis/SpinComponent.hpp"
#include "Basis/TransformationMatrix.hpp"
#include "Mathematical/Representation/QCMatrix.hpp"
#include "Mathematical/ScalarFunction.hpp"
#include "Processing/RDM/OneRDM.hpp"
#include "Processing/RDM/TwoRDM.hpp"
#include "Utilities/type_traits.hpp"

#include <array>


namespace GQCP {


/**
 *  A class that represents an unrestricted second-quantized one-electron operator: it holds the matrix representation of its parameters for both spin components, which are (usually) integrals over first-quantized operators
 *
 *  @tparam _Scalar             the scalar type, i.e. the scalar representation of one of the parameters
 *  @tparam _Components         the number of components of the second-quantized operator
 */
template <typename _Scalar, size_t _Components>
class USQOneElectronOperator {
public:
    using Scalar = _Scalar;
    static constexpr auto Components = _Components;


private:
    std::array<QCMatrix<Scalar>, Components> fs_alpha;  // all the matrix representations of the alpha spin component (hence the s) of the parameters (integrals) of the different components of this second-quantized operator
    std::array<QCMatrix<Scalar>, Components> fs_beta;  // all the matrix representations of the beta spin component (hence the s) of the parameters (integrals) of the different components of this second-quantized operator


public:

    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param fs_alpha    all the matrix representations of the alpha spin component (hence the s) of the parameters (integrals) of the different components of this second-quantized operator
     *  @param fs_beta     all the matrix representations of the beta spin component (hence the s) of the parameters (integrals) of the different components of this second-quantized operator
     */
    USQOneElectronOperator(const std::array<QCMatrix<Scalar>, Components>& fs_alpha, const std::array<QCMatrix<Scalar>, Components>& fs_beta) :
        fs_alpha (fs_alpha),
        fs_beta (fs_beta)
    {
        // Check if the given matrix representations have the same dimensions.
        const auto dimension_of_first_alpha = this->fs_alpha[0].dimension();
        const auto dimension_of_first_beta = this->fs_beta[0].dimension();

        for (size_t i = 1; i < Components; i++) {

            const auto dimension_of_ith_alpha = this->fs_alpha[i].dimension();
            const auto dimension_of_ith_beta = this->fs_beta[i].dimension();

            if ((dimension_of_first_alpha != dimension_of_ith_alpha) || (dimension_of_first_beta != dimension_of_ith_beta)) {
                throw std::invalid_argument("USQOneElectronOperator(const std::array<QCMatrix<Scalar>, Components>&, const std::array<QCMatrix<Scalar>, components>&): The given matrix representations do not have the same dimensions for either the alpha or beta component.");
            }
        }
    }


    /**
     *  A constructor for ScalarUSQOneElectronOperators that doesn't require the argument to be a vector of just one element.
     * 
     *  @param f_alpha            the matrix representation of the integrals of this scalar second-quantized operator, for the alpha spin component
     *  @param f_beta             the matrix representation of the integrals of this scalar second-quantized operator, for the beta spin component
     * 
     *  @note This constructor is only available for ScalarUSQOneElectronOperators (for the std::enable_if, see https://stackoverflow.com/a/17842695/7930415)
     */
    template <size_t Z = Components>
    USQOneElectronOperator(const QCMatrix<Scalar>& f_alpha, const QCMatrix<Scalar>& f_beta, typename std::enable_if<Z == 1>::type* = 0) :
        USQOneElectronOperator(std::array<QCMatrix<Scalar>, 1>{f_alpha}, std::array<QCMatrix<Scalar>, 1>{f_beta})
    {}


    /**
     *  Construct a one-electron operator with parameters that are zero, for both the alpha and beta spin components
     * 
     *  @param dim_alpha          the dimension of the alpha matrix representation of the parameters, i.e. the number of orbitals/sites
     *  @param dim_beta           the dimension of the beta matrix representation of the parameters, i.e. the number of orbitals/sites
     * 
     */
    USQOneElectronOperator(const size_t dim_alpha, const size_t dim_beta) {
        for (size_t i = 0; i < Components; i++) {
            this->fs_alpha[i] = QCMatrix<Scalar>::Zero(dim_alpha, dim_alpha);
            this->fs_beta[i] = QCMatrix<Scalar>::Zero(dim_beta, dim_beta);
        }
    }


    /**
     *  Construct a one-electron operator with parameters that are zero, for both the alpha and beta spin components
     * 
     *  @param dim          the dimension of the matrix representation of the parameters, i.e. the number of orbitals/sites. This dimension is the same for the alpha and beta component.
     */
    USQOneElectronOperator(const size_t dim) :
        USQOneElectronOperator(dim, dim)
    {}


    /**
     *  Default constructor: construct a one-electron operator with parameters that are zero, for both spin components.
     */
    USQOneElectronOperator() :
        USQOneElectronOperator(0)  // dimensions of the representations are zero
    {}


    /*
     *  OPERATORS
     */

    /**
     *  @param i            the index
     * 
     *  @return the i-th component (i.e. the i-th alpha and i-th beta component) of this operator
     */
    USQOneElectronOperator<Scalar, 1> operator[](const size_t i) const {

        if (i >= Components) {
            throw std::invalid_argument("USQOneElectronOperator::operator[](const size_t): The given index is out of bounds.");
        }

        return USQOneElectronOperator<Scalar, 1> {this->fs_alpha[i], this->fs_beta[i]};
    }


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return read-only matrix representations of all the alpha parameters (integrals) of the different components of this second-quantized operator
     */
    const std::array<QCMatrix<Scalar>, Components>& allAlphaParameters() const {
        return this->fs_alpha;
    }


    /**
     *  @return writable matrix representations of all the alpha parameters (integrals) of the different components of this second-quantized operator
     */
    std::array<QCMatrix<Scalar>, Components>& allAlphaParameters() {
        return this->fs_alpha;
    }


    /**
     *  @return read-only matrix representations of all the beta parameters (integrals) of the different components of this second-quantized operator
     */
    const std::array<QCMatrix<Scalar>, Components>& allBetaParameters() const {
        return this->fs_beta;
    }


    /**
     *  @return writable matrix representations of all the beta parameters (integrals) of the different components of this second-quantized operator
     */
    std::array<QCMatrix<Scalar>, Components>& allBetaParameters() {
        return this->fs_beta;
    }


    /**
     *  @return the dimension of the alpha matrices
     */
    size_t alphaDimension() const { return this->fs_alpha[0].dimension(); }


    /**
     *  @param i            the index of the component
     * 
     *  @return a read-only the matrix representation of the parameters (integrals) of one of the the different alpha components of this second-quantized operator
     */
    const QCMatrix<Scalar>& alphaParameters(const size_t i = 0) const {
        return this->fs_alpha[i];
    }


    /**
     *  @param i            the index of the component
     * 
     *  @return a read-only the matrix representation of the parameters (integrals) of one of the the different alpha components of this second-quantized operator
     */
    QCMatrix<Scalar>& alphaParameters(const size_t i = 0) {
        return this->fs_alpha[i];
    }


    /**
     *  @return the dimension of the beta matrices
     */
    size_t betaDimension() const { return this->fs_beta[0].dimension(); }


    /**
     *  @param i            the index of the component
     * 
     *  @return a read-only the matrix representation of the parameters (integrals) of one of the the different beta components of this second-quantized operator
     */
    const QCMatrix<Scalar>& betaParameters(const size_t i = 0) const {
        return this->fs_beta[i];
    }


    /**
     *  @param i            the index of the component
     * 
     *  @return a read-only the matrix representation of the parameters (integrals) of one of the the different beta components of this second-quantized operator
     */
    QCMatrix<Scalar>& betaParameters(const size_t i = 0) {
        return this->fs_beta[i];
    }


    /**
     *  @param D_alpha                the alpha 1-RDM that represents the wave function
     *  @param D_beta                 the beta 1-RDM that represents the wave function
     *
     *  @return the expectation values of all components of the one-electron operator
     */
    Vector<Scalar, Components> calculateExpectationValue(const OneRDM<Scalar>& D_alpha, const OneRDM<Scalar>& D_beta) const {

        if (this->fs_alpha[0].dimension() != D_alpha.dimension() || this->fs_beta[0].dimension() != D_beta.dimension()) {
            throw std::invalid_argument("USQOneElectronOperator::calculateExpectationValue(const OneRDM<Scalar>, OneRDM<Scalar>): The given 1-RDM is not compatible with the one-electron operator.");
        }

        std::array<Scalar, Components> expectation_values {} ; // zero initialization

        for (size_t i = 0; i < Components; i++) {
            expectation_values[i] = (this->alphaParameters(i) * D_alpha).trace() + (this->betaParameters(i) * D_beta).trace();
        }

        return Eigen::Map<Eigen::Matrix<Scalar, Components, 1>>(expectation_values.data());  // convert std::array to Vector
    }


    /**
     *  @return the sum of the alpha and beta dimensions
     */
    size_t dimension() const {
        return this->alphaDimension() + this->betaDimension();
    }


    /**
     *  In-place rotate the operator to another basis. The same matrix is used for the alpha and beta components.   
     * 
     *  @param U                            the (unitary) rotation matrix
     */
    void rotate(const TransformationMatrix<Scalar>& U) {

        // Transform the matrix representations of the components
        for (auto& f_alpha : this->allAlphaParameters()) {
            f_alpha.basisRotateInPlace(U);
        }
        for (auto& f_beta : this->allBetaParameters()) {
            f_beta.basisRotateInPlace(U);
        }
    }


    /**
     *  In-place rotate the operator using a unitary Jacobi rotation matrix constructed from the Jacobi rotation parameters. The same matrix is used to rotate the alpha and beta components.
     * 
     *  @param jacobi_rotation_parameters       the Jacobi rotation parameters (p, q, angle) that are used to specify a Jacobi rotation: we use the (cos, sin, -sin, cos) definition for the Jacobi rotation matrix
     */
    void rotate(const JacobiRotationParameters& jacobi_rotation_parameters) {

        // Transform the matrix representations of the components
        for (auto& f_alpha : this->allAlphaParameters()) {
            f_alpha.basisRotateInPlace(jacobi_rotation_parameters);
        }
        for (auto& f_beta : this->allBetaParameters()) {
            f_beta.basisRotateInPlace(jacobi_rotation_parameters);
        }
    }


    /**
     *  In-place transform the operator to another basis. The same transformation is applied on the alpha and beta components.
     * 
     *  @param T                            the transformation matrix
     */
    void transform(const TransformationMatrix<Scalar>& T) {

        // Transform the matrix representations of the components
        for (auto& f_alpha : this->allAlphaParameters()) {
            f_alpha.basisTransformInPlace(T);
        }
        for (auto& f_beta : this->allBetaParameters()) {
            f_beta.basisTransformInPlace(T);
        }
    }
};



/*
 *  CONVENIENCE ALIASES
 */
template <typename Scalar>
using ScalarUSQOneElectronOperator = USQOneElectronOperator<Scalar, 1>;

template <typename Scalar>
using VectorUSQOneElectronOperator = USQOneElectronOperator<Scalar, 3>;


/*
 *  OPERATORS
 */

/**
 *  Add two one-electron operators by adding their parameters. The two alphas are added together and the two betas are added together. 
 * 
 *  @tparam LHSScalar           the scalar type of the left-hand side
 *  @tparam RHSScalar           the scalar type of the right-hand side
 *  @tparam Components          the number of components of the one-electron operators
 * 
 *  @param lhs                  the left-hand side
 *  @param rhs                  the right-hand side
 */
template <typename LHSScalar, typename RHSScalar, size_t Components>
auto operator+(const USQOneElectronOperator<LHSScalar, Components>& lhs, const USQOneElectronOperator<RHSScalar, Components>& rhs) -> USQOneElectronOperator<sum_t<LHSScalar, RHSScalar>, Components> {

    using ResultScalar = sum_t<LHSScalar, RHSScalar>;

    auto F_sum_alpha = lhs.allAlphaParameters();
    auto F_sum_beta = lhs.allBetaParameters();
    for (size_t i = 0; i < Components; i++) {
        F_sum_alpha[i] += rhs.alphaParameters(i);
        F_sum_beta[i] += rhs.betaParameters(i);
    }

    return USQOneElectronOperator<ResultScalar, Components>(F_sum_alpha, F_sum_beta);
}


/**
 *  Multiply a one-electron operator with a scalar. The alpha and beta components are multiplied with the same scalar.
 * 
 *  @tparam Scalar              the scalar type of the scalar
 *  @tparam OperatorScalar      the scalar type of the operator
 * 
 *  @tparam scalar              the scalar of the scalar multiplication
 *  @tparam op                  the one-electron operator
 */
template <typename Scalar, typename OperatorScalar, size_t Components>
auto operator*(const Scalar& scalar, const USQOneElectronOperator<OperatorScalar, Components>& op) -> USQOneElectronOperator<product_t<Scalar, OperatorScalar>, Components> {

    using ResultScalar = product_t<Scalar, OperatorScalar>;

    auto fs_alpha = op.allAlphaParameters();
    auto fs_beta = op.allBetaParameters();
    for (auto& f_alpha : fs_alpha) {
        f_alpha *= scalar;
    }
    for (auto& f_beta : fs_beta) {
        f_beta *= scalar;
    }

    return USQOneElectronOperator<ResultScalar, Components>(fs_alpha, fs_beta);
}


/**
 *  Negate a one-electron operator
 * 
 *  @tparam Scalar              the scalar type of the operator
 *  @tparam Components          the number of components of the one-electron operator
 * 
 *  @param op                   the operator
 */
template <typename Scalar, size_t Components>
USQOneElectronOperator<Scalar, Components> operator-(const USQOneElectronOperator<Scalar, Components>& op) {

    return (-1.0) * op;  // negation is scalar multiplication with (-1.0)
}


/**
 *  Subtract two one-electron operators by subtracting their parameters. The alphas are subtracted from the alphas and the betas from the betas.
 * 
 *  @tparam LHSScalar           the scalar type of the left-hand side
 *  @tparam RHSScalar           the scalar type of the right-hand side
 *  @tparam Components          the number of components of the one-electron operators
 * 
 *  @param lhs                  the left-hand side
 *  @param rhs                  the right-hand side
 */
template <typename LHSScalar, typename RHSScalar, size_t Components>
auto operator-(const USQOneElectronOperator<LHSScalar, Components>& lhs, const USQOneElectronOperator<RHSScalar, Components>& rhs) -> USQOneElectronOperator<sum_t<LHSScalar, RHSScalar>, Components> {

    using ResultScalar = sum_t<LHSScalar, RHSScalar>;

    auto F_min_alpha = lhs.allAlphaParameters();
    auto F_min_beta = lhs.allBetaParameters();
    for (size_t i = 0; i < Components; i++) {
        F_min_alpha[i] -= rhs.alphaParameters(i);
        F_min_beta[i] -= rhs.betaParameters(i);
    }

    return USQOneElectronOperator<ResultScalar, Components>(F_min_alpha, F_min_beta);
}


} // namespace GQCP
