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
#include "Basis/TransformationMatrix.hpp"
#include "Mathematical/Representation/QCRankFourTensor.hpp"
#include "Operator/SecondQuantized/USQOneElectronOperator.hpp"
#include "Processing/RDM/OneRDM.hpp"
#include "Processing/RDM/TwoRDM.hpp"
#include "Utilities/miscellaneous.hpp"

#include <array>


namespace GQCP {


/**
A class that represents an 'unrestricted second-quantized two-electron operator' suitable for the projection of the non-relativistic Hamiltonian onto an unrestricted spinor basis. It holds the tensor representation of its parameters for both spin components and both mixed spin components, which are (usually) integrals over first-quantized operators
 * 
 *  @tparam _Scalar             the scalar type, i.e. the scalar representation of one of the parameters
 *  @tparam _Components         the number of components of the second-quantized operator
 */
template <typename _Scalar, size_t _Components>
class USQTwoElectronOperator {
public:

    using Scalar = _Scalar;
    static constexpr auto Components = _Components;


private:
    std::array<QCRankFourTensor<Scalar>, Components> gs_alphaAlpha;  // all the tensor representations (hence the 's') of the alpha-alpha parameters (integrals) of the different components of this second-quantized operator
    std::array<QCRankFourTensor<Scalar>, Components> gs_alphaBeta;  // all the tensor representations (hence the 's') of the alpha-beta parameters (integrals) of the different components of this second-quantized operator
    std::array<QCRankFourTensor<Scalar>, Components> gs_betaAlpha;  // all the tensor representations (hence the 's') of the beta-alpha parameters (integrals) of the different components of this second-quantized operator
    std::array<QCRankFourTensor<Scalar>, Components> gs_betaBeta;  // all the tensor representations (hence the 's') of the beta-beta parameters (integrals) of the different components of this second-quantized operator

public:

    /*
     *  CONSTRUCTORS
     */
    /**
     *  @param gs_alphaAlpha            all the tensor representations (hence the 's') of the alpha-alpha parameters (integrals) of the different components of this second-quantized operator
     *  @param gs_alphaBeta             all the tensor representations (hence the 's') of the alpha-beta parameters (integrals) of the different components of this second-quantized operator
     *  @param gs_betaAlpha             all the tensor representations (hence the 's') of the beta-alpha parameters (integrals) of the different components of this second-quantized operator
     *  @param gs_betaBeta              all the tensor representations (hence the 's') of the beta-beta parameters (integrals) of the different components of this second-quantized operator
     * 
     */
    USQTwoElectronOperator(const std::array<QCRankFourTensor<Scalar>, Components>& gs_alphaAlpha, const std::array<QCRankFourTensor<Scalar>, Components>& gs_alphaBeta, const std::array<QCRankFourTensor<Scalar>, Components>& gs_betaAlpha, const std::array<QCRankFourTensor<Scalar>, Components>& gs_betaBeta) : 
        gs_alphaAlpha (gs_alphaAlpha),
        gs_alphaBeta (gs_alphaBeta),
        gs_betaAlpha (gs_betaAlpha),
        gs_betaBeta (gs_betaBeta)
    {
        // Check if the given tensor representations have the same dimensions
        const auto dimension_of_first_alphaAlpha = this->gs_alphaAlpha[0].dimension();
        const auto dimension_of_first_alphaBeta = this-> gs_alphaBeta[0].dimension();
        const auto dimension_of_first_betaAlpha = this->gs_betaAlpha[0].dimension();
        const auto dimension_of_first_betaBeta = this-> gs_betaBeta[0].dimension();

        for (size_t i = 1; i < Components; i++) {

            const auto dimension_of_ith_alphaAlpha = this->gs_alphaAlpha[i].dimension();
            const auto dimension_of_ith_alphaBeta = this->gs_alphaBeta[i].dimension();
            const auto dimension_of_ith_betaAlpha = this->gs_betaAlpha[i].dimension();
            const auto dimension_of_ith_betaBeta = this->gs_betaBeta[i].dimension();
            if ((dimension_of_first_alphaAlpha != dimension_of_ith_alphaAlpha) || (dimension_of_first_alphaBeta != dimension_of_ith_alphaBeta) || (dimension_of_first_betaAlpha != dimension_of_ith_betaAlpha) || (dimension_of_first_betaBeta != dimension_of_ith_betaBeta)) {
                throw std::invalid_argument("SQTwoElectronOperator(const std::array<QCMatrix<Scalar>, Components>&): The given matrix representations do not have the same dimensions for either the alpha or beta component.");
            }
        }
    }


    /**
     *  A constructor for ScalarUSQTwoElectronOperators that doesn't require the argument to be a vector of just one element.
     * 
     *  @param g_alphaAlpha           the tensor representation of the alpha-alpha integrals of this scalar second-quantized operator
     *  @param g_alphaBeta            the tensor representation of the alpha-beta integrals of this scalar second-quantized operator
     *  @param g_betaAlpha            the tensor representation of the beta-alpha integrals of this scalar second-quantized operator
     *  @param g_betaBEta             the tensor representation of the beta-beta integrals of this scalar second-quantized operator
     * 
     *  @note This constructor is only available for ScalarSQTwoElectronOperators (for the std::enable_if, see https://stackoverflow.com/a/17842695/7930415)
     */
    template <size_t Z = Components>
    USQTwoElectronOperator(const QCRankFourTensor<Scalar>& g_alphaAlpha, const QCRankFourTensor<Scalar>& g_alphaBeta, const QCRankFourTensor<Scalar>& g_betaAlpha, const QCRankFourTensor<Scalar>& g_betaBeta, typename std::enable_if<Z == 1>::type* = 0) :
        USQTwoElectronOperator(std::array<QCRankFourTensor<Scalar>, 1>{g_alphaAlpha}, std::array<QCRankFourTensor<Scalar>, 1>{g_alphaBeta}, std::array<QCRankFourTensor<Scalar>, 1>{g_betaAlpha}, std::array<QCRankFourTensor<Scalar>, 1>{g_betaBeta})
    {}


    /**
     *  Construct an unrestricted two-electron operator with zero parameters, dimensions of alpha and beta component are the same.
     * 
     *  @param dim        the dimension of the matrix representation of the alpha and beta parameters, i.e. the number of orbitals/sites
     * 
     */
    USQTwoElectronOperator(const size_t dim) {
        for (size_t i = 0; i < Components; i++) {
            this->gs_alphaAlpha[i] = QCRankFourTensor<Scalar>(dim);
            this->gs_alphaBeta[i] = QCRankFourTensor<Scalar>(dim);
            this->gs_betaAlpha[i] = QCRankFourTensor<Scalar>(dim);
            this->gs_betaBeta[i] = QCRankFourTensor<Scalar>(dim);
            this->gs_alphaAlpha[i].setZero();
            this->gs_alphaBeta[i].setZero();
            this->gs_betaAlpha[i].setZero();
            this->gs_betaBeta[i].setZero();
        }
    }

    
    /**
     *  Default constructor: construct an unrestricted two-electron operator with parameters that are zero
     */
    USQTwoElectronOperator() :
        USQTwoElectronOperator(0)  // dimensions of the representations are zero
    {}


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return read-only tensor representations of all the alpha-alpha parameters (integrals) of the different components of this second-quantized operator
     */
    const std::array<QCRankFourTensor<Scalar>, Components>& allAlphaAlphaParameters() const {
        return this->gs_alphaAlpha;
    }

    /**
     *  @return writable tensor representations of all the alpha-alpha parameters (integrals) of the different components of this second-quantized operator
     */
    std::array<QCRankFourTensor<Scalar>, Components>& allAlphaAlphaParameters() {
        return this->gs_alphaAlpha;
    }


    /**
     *  @return read-only tensor representations of all the alpha-beta parameters (integrals) of the different components of this second-quantized operator
     */
    const std::array<QCRankFourTensor<Scalar>, Components>& allAlphaBetaParameters() const {
        return this->gs_alphaBeta;
    }

    /**
     *  @return writable tensor representations of all the alpha-beta parameters (integrals) of the different components of this second-quantized operator
     */
    std::array<QCRankFourTensor<Scalar>, Components>& allAlphaBetaParameters() {
        return this->gs_alphaBeta;
    }


    /**
     *  @return read-only tensor representations of all the beta-alpha parameters (integrals) of the different components of this second-quantized operator
     */
    const std::array<QCRankFourTensor<Scalar>, Components>& allBetaAlphaParameters() const {
        return this->gs_betaAlpha;
    }


    /**
     *  @return writable tensor representations of all the beta-alpha parameters (integrals) of the different components of this second-quantized operator
     */
    std::array<QCRankFourTensor<Scalar>, Components>& allBetaAlphaParameters() {
        return this->gs_betaAlpha;
    }


    /**
     *  @return read-only tensor representations of all the beta-beta parameters (integrals) of the different components of this second-quantized operator
     */
    const std::array<QCRankFourTensor<Scalar>, Components>& allBetaBetaParameters() const {
        return this->gs_betaBeta;
    }


    /**
     *  @return writable tensor representations of all the beta-beta parameters (integrals) of the different components of this second-quantized operator
     */
    std::array<QCRankFourTensor<Scalar>, Components>& allBetaBetaParameters() {
        return this->gs_betaBeta;
    }


    /** partition
      *  @return the dimension of the alpha-alpha components
     */
    size_t alphaAlphaDimension() const { return this->gs_alphaAlpha[0].dimension(); }


    /**
     *  @return the dimension of the alpha-beta components
     */
    size_t alphaBetaDimension() const { return this->gs_alphaBeta[0].dimension(); }


    /**
     *  @param i            the index of the component
     * 
     *  @return a read-only tensor representation of the alpha-alpha parameters (integrals) of one of the the different components of this second-quantized operator
     */
    const QCRankFourTensor<Scalar>& alphaAlphaParameters(const size_t i = 0) const {
        return this->gs_alphaAlpha[i];
    }


    /**
     *  @param i            the index of the component
     * 
     *  @return a writable tensor representation of the alpha-alpha parameters (integrals) of one of the the different components of this second-quantized operator
     */
    QCRankFourTensor<Scalar>& alphaAlphaParameters(const size_t i = 0) {
        return this->gs_alphaAlpha[i];
    }


    /**
     *  @param i            the index of the component
     * 
     *  @return a read-only tensor representation of the alpha-beta parameters (integrals) of one of the the different components of this second-quantized operator
     */
    const QCRankFourTensor<Scalar>& alphaBetaParameters(const size_t i = 0) const {
        return this->gs_alphaBeta[i];
    }


    /**
     *  @param i            the index of the component
     * 
     *  @return a writable tensor representation of the alpha-beta parameters (integrals) of one of the the different components of this second-quantized operator
     */
    QCRankFourTensor<Scalar>& alphaBetaParameters(const size_t i = 0) {
        return this->gs_alphaBeta[i];
    }


    /**
     *  @return the dimension of the beta-alpha components
     */
    size_t betaAlphaDimension() const { return this->gs_betaAlpha[0].dimension(); }


    /**
     *  @param i            the index of the component
     * 
     *  @return a read-only tensor representation of the beta-alpha parameters (integrals) of one of the the different components of this second-quantized operator
     */
    const QCRankFourTensor<Scalar>& betaAlphaParameters(const size_t i = 0) const {
        return this->gs_betaAlpha[i];
    }


    /**
     *  @param i            the index of the component
     * 
     *  @return a writable tensor representation of the beta-alpha parameters (integrals) of one of the the different components of this second-quantized operator
     */
    QCRankFourTensor<Scalar>& betaAlphaParameters(const size_t i = 0) {
        return this->gs_betaAlpha[i];
    }


    /**
     *  @return the dimension of the beta-beta components
     */
    size_t betaBetaDimension() const { return this->gs_betaBeta[0].dimension(); }


    /**
     *  @param i            the index of the component
     * 
     *  @return a read-only tensor representation of the beta-beta parameters (integrals) of one of the the different components of this second-quantized operator
     */
    const QCRankFourTensor<Scalar>& betaBetaParameters(const size_t i = 0) const {
        return this->gs_betaBeta[i];
    }


    /**
     *  @param i            the index of the component
     * 
     *  @return a writable tensor representation of the beta-beta parameters (integrals) of one of the the different components of this second-quantized operator
     */
    QCRankFourTensor<Scalar>& betaBetaParameters(const size_t i = 0) {
        return this->gs_betaBeta[i];
    }

    
    /**
     *  @param d            the 2-RDM that represents the wave function
     *
     *  @return the expectation values of all the components of the two-electron operator, with the given 2-RDM: this includes the prefactor 1/2
     */
    Vector<Scalar, Components> calculateExpectationValue(const TwoRDM<Scalar>& d_alphaAlpha, const TwoRDM<Scalar>& d_alphaBeta, const TwoRDM<Scalar>& d_betaAlpha, const TwoRDM<Scalar>& d_betaBeta) const {

        if ((this->alphaAlphaDimension() != d_alphaAlpha.dimension()) || (this->alphaBetaDimension() != d_alphaBeta.dimension()) || (this->betaAlphaDimension() != d_betaAlpha.dimension()) || (this->betaBetaDimension() != d_betaBeta.dimension())) {
            throw std::invalid_argument("SQTwoElectronOperator::calculateExpectationValue(const TwoRDM<double>&): The given 2-RDM is not compatible with the respective component of the two-electron operator.");
        }


        std::array<Scalar, Components> expectation_values {};  // zero initialization of alphha component
        for (size_t i = 0; i < Components; i++) {

            // Specify the contractions for the relevant contraction of the two-electron integrals and the 2-RDM
            //      0.5 g(p q r s) d(p q r s)
            Eigen::array<Eigen::IndexPair<int>, 4> contractions = {Eigen::IndexPair<int>(0,0), Eigen::IndexPair<int>(1,1), Eigen::IndexPair<int>(2,2), Eigen::IndexPair<int>(3,3)};
            //      Perform the contraction
            Eigen::Tensor<Scalar, 0> contraction_alphaAlpha = 0.5 * this->alphaAlphaParameters(i).contract(d_alphaAlpha.Eigen(), contractions);
            Eigen::Tensor<Scalar, 0> contraction_alphaBeta = 0.5 * this->alphaBetaParameters(i).contract(d_alphaBeta.Eigen(), contractions);
            Eigen::Tensor<Scalar, 0> contraction_betaAlpha = 0.5 * this->betaAlphaParameters(i).contract(d_betaAlpha.Eigen(), contractions);
            Eigen::Tensor<Scalar, 0> contraction_betaBeta = 0.5 * this->betaBetaParameters(i).contract(d_betaBeta.Eigen(), contractions);

            // As the contraction is a scalar (a tensor of rank 0), we should access by (0).
            expectation_values[i] = contraction_alphaAlpha(0) + contraction_alphaBeta(0) + contraction_betaAlpha(0) + contraction_betaBeta(0);
        }

        return Eigen::Map<Eigen::Matrix<Scalar, Components, 1>>(expectation_values.data());  // convert std::array to Vector
    }


    /**
     *  In-place rotate the operator to another basis. The alpha-alpha, alpha-beta, beta-alpha and beta-beta components are transformed in the same way.
     * 
     *  @param U                            the (unitary) rotation matrix
     */
    void rotate(const TransformationMatrix<Scalar>& U) {

        // Transform the matrix representations of the components
        for (auto& g_alphaAlpha : this->allAlphaAlphaParameters()) {
            g_alphaAlpha.basisRotateInPlace(U);
        }
        for (auto& g_alphaBeta : this->allAlphaBetaParameters()) {
            g_alphaBeta.basisRotateInPlace(U);
        }
        for (auto& g_betaAlpha : this->allBetaAlphaParameters()) {
            g_betaAlpha.basisRotateInPlace(U);
        }
        for (auto& g_betaBeta : this->allBetaBetaParameters()) {
            g_betaBeta.basisRotateInPlace(U);
        }
    }


    /**
     *  In-place rotate the operator using a unitary Jacobi rotation matrix constructed from the Jacobi rotation parameters. The alpha-alpha, alpha-beta, beta-alpha and beta-beta components are transformed in the same way.
     * 
     *  @param jacobi_rotation_parameters       the Jacobi rotation parameters (p, q, angle) that are used to specify a Jacobi rotation: we use the (cos, sin, -sin, cos) definition for the Jacobi rotation matrix
     */
    void rotate(const JacobiRotationParameters& jacobi_rotation_parameters) {

        // Transform the matrix representations of the components
        for (auto& g_alphaAlpha : this->allAlphaAlphaParameters()) {
            g_alphaAlpha.basisRotateInPlace(jacobi_rotation_parameters);
        }
        for (auto& g_alphaBeta : this->allAlphaBetaParameters()) {
            g_alphaBeta.basisRotateInPlace(jacobi_rotation_parameters);
        }
        for (auto& g_betaAlpha : this->allBetaAlphaParameters()) {
            g_betaAlpha.basisRotateInPlace(jacobi_rotation_parameters);
        }
        for (auto& g_betaBeta : this->allBetaBetaParameters()) {
            g_betaBeta.basisRotateInPlace(jacobi_rotation_parameters);
        }
    }


    /**
     *  In-place transform the operator to another basis. The alpha-alpha, alpha-beta, beta-alpha and beta-beta components are transformed in the same way.
     * 
     *  @param T                            the transformation matrix
     */
    void transform(const TransformationMatrix<Scalar>& T) {

        // Transform the matrix representations of the components
        for (auto& g_alphaAlpha : this->allAlphaAlphaParameters()) {
            g_alphaAlpha.basisTransformInPlace(T);
        }
        for (auto& g_alphaBeta : this->allAlphaBetaParameters()) {
            g_alphaBeta.basisTransformInPlace(T);
        }
        for (auto& g_betaAlpha : this->allBetaAlphaParameters()) {
            g_betaAlpha.basisTransformInPlace(T);
        }
        for (auto& g_betaBeta : this->allBetaBetaParameters()) {
            g_betaBeta.basisTransformInPlace(T);
        }
    }
};



/*
 *  CONVENIENCE ALIASES
 */
template <typename Scalar>
using ScalarUSQTwoElectronOperator = USQTwoElectronOperator<Scalar, 1>;


/*
 *  OPERATORS
 */

/**
 *  Add two two-electron operators by adding their parameters. The two alpha-alpha components are added together and the respective alpha-beta, beta-alpha and beta-beta components are added analogously.
 * 
 *  @tparam LHSScalar           the scalar type of the left-hand side
 *  @tparam RHSScalar           the scalar type of the right-hand side
 *  @tparam Components          the number of components of the two-electron operators
 * 
 *  @param lhs                  the left-hand side
 *  @param rhs                  the right-hand side
 */
template <typename LHSScalar, typename RHSScalar, size_t Components>
auto operator+(const USQTwoElectronOperator<LHSScalar, Components>& lhs, const USQTwoElectronOperator<RHSScalar, Components>& rhs) -> USQTwoElectronOperator<sum_t<LHSScalar, RHSScalar>, Components> {

    using ResultScalar = sum_t<LHSScalar, RHSScalar>;

    auto G_sum_alphaAlpha = lhs.allAlphaAlphaParameters();
    auto G_sum_alphaBeta = lhs.allAlphaBetaParameters();
    auto G_sum_betaAlpha = lhs.allBetaAlphaParameters();
    auto G_sum_betaBeta = lhs.allBetaBetaParameters();
    for (size_t i = 0; i < Components; i++) {
        G_sum_alphaAlpha[i] += rhs.alphaAlphaParameters(i);
        G_sum_alphaBeta[i] += rhs.alphaBetaParameters(i);
        G_sum_betaAlpha[i] += rhs.betaAlphaParameters(i);
        G_sum_betaBeta[i] += rhs.betaBetaParameters(i);
    }

    return USQTwoElectronOperator<ResultScalar, Components>(G_sum_alphaAlpha, G_sum_alphaBeta, G_sum_betaAlpha, G_sum_betaBeta);
}


/**
 *  Multiply a two-electron operator with a scalar. The alpha-alpha, alpha-beta, beta-alpha and beta-beta components are multiplied with the same scalar.
 * 
 *  @tparam Scalar              the scalar type of the scalar
 *  @tparam OperatorScalar      the scalar type of the operator
 * 
 *  @tparam scalar              the scalar of the scalar multiplication
 *  @tparam op                  the two-electron operator
 */
template <typename Scalar, typename OperatorScalar, size_t Components>
auto operator*(const Scalar& scalar, const USQTwoElectronOperator<OperatorScalar, Components>& op) -> USQTwoElectronOperator<product_t<Scalar, OperatorScalar>, Components> {

    using ResultScalar = product_t<Scalar, OperatorScalar>;

    auto G_alphaAlpha = op.allAlphaAlphaParameters();
    auto G_alphaBeta = op.allAlphaBetaParameters();
    auto G_betaAlpha = op.allBetaAlphaParameters();
    auto G_betaBeta = op.allBetaBetaParameters();
    for (size_t i = 0; i < Components; i++) {
        G_alphaAlpha[i] = scalar * G_alphaAlpha[i];
        G_alphaBeta[i] = scalar * G_alphaBeta[i];
        G_betaAlpha[i] = scalar * G_betaAlpha[i];
        G_betaBeta[i] = scalar * G_betaBeta[i];
    }

    return USQTwoElectronOperator<ResultScalar, Components>(G_alphaAlpha, G_alphaBeta, G_betaAlpha, G_betaBeta);
}


/**
 *  Negate a two-electron operator
 * 
 *  @tparam Scalar              the scalar type of the operator
 *  @tparam Components          the number of components of the one-electron operator
 * 
 *  @param op                   the operator
 */
template <typename Scalar, size_t Components>
USQTwoElectronOperator<Scalar, Components> operator-(const USQTwoElectronOperator<Scalar, Components>& op) {

    return (-1.0) * op;  // negation is scalar multiplication with (-1.0)
}


/**
 *  Subtract two two-electron operators by subtracting their respective alpha-alpha, alpha-beta, beta-alpha and beta-beta parameters
 * 
 *  @tparam LHSScalar           the scalar type of the left-hand side
 *  @tparam RHSScalar           the scalar type of the right-hand side
 *  @tparam Components          the number of components of the two-electron operators
 * 
 *  @param lhs                  the left-hand side
 *  @param rhs                  the right-hand side
 */
template <typename LHSScalar, typename RHSScalar, size_t Components>
auto operator-(const USQTwoElectronOperator<LHSScalar, Components>& lhs, const USQTwoElectronOperator<RHSScalar, Components>& rhs) -> USQTwoElectronOperator<sum_t<LHSScalar, RHSScalar>, Components> {

    return lhs + (-rhs);
}


} // namespace GQCP
