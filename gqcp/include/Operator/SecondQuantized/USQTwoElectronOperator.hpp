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
    std::array<QCRankFourTensor<Scalar>, Components> gs_aa;  // all the tensor representations (hence the 's') of the alpha-alpha parameters (integrals) of the different components of this second-quantized operator
    std::array<QCRankFourTensor<Scalar>, Components> gs_ab;  // all the tensor representations (hence the 's') of the alpha-beta parameters (integrals) of the different components of this second-quantized operator
    std::array<QCRankFourTensor<Scalar>, Components> gs_ba;  // all the tensor representations (hence the 's') of the beta-alpha parameters (integrals) of the different components of this second-quantized operator
    std::array<QCRankFourTensor<Scalar>, Components> gs_bb;  // all the tensor representations (hence the 's') of the beta-beta parameters (integrals) of the different components of this second-quantized operator

public:

    /*
     *  CONSTRUCTORS
     */
    /**
     *  @param gs_aa            all the tensor representations (hence the 's') of the alpha-alpha parameters (integrals) of the different components of this second-quantized operator
     *  @param gs_ab            all the tensor representations (hence the 's') of the alpha-beta parameters (integrals) of the different components of this second-quantized operator
     *  @param gs_ba            all the tensor representations (hence the 's') of the beta-alpha parameters (integrals) of the different components of this second-quantized operator
     *  @param gs_bb            all the tensor representations (hence the 's') of the beta-beta parameters (integrals) of the different components of this second-quantized operator
     * 
     */
    USQTwoElectronOperator(const std::array<QCRankFourTensor<Scalar>, Components>& gs_aa, const std::array<QCRankFourTensor<Scalar>, Components>& gs_ab, const std::array<QCRankFourTensor<Scalar>, Components>& gs_ba, const std::array<QCRankFourTensor<Scalar>, Components>& gs_bb) : 
        gs_aa (gs_aa),
        gs_ab (gs_ab),
        gs_ba (gs_ba),
        gs_bb (gs_bb)
    {
        // Check if the given tensor representations have the same dimensions, for each spin part.
        const auto dimension_of_first_aa = this->gs_aa[0].dimension();
        const auto dimension_of_first_ab = this-> gs_ab[0].dimension();
        const auto dimension_of_first_ba = this->gs_ba[0].dimension();
        const auto dimension_of_first_bb = this-> gs_bb[0].dimension();

        for (size_t i = 1; i < Components; i++) {

            const auto dimension_of_ith_aa = this->gs_aa[i].dimension();
            const auto dimension_of_ith_ab = this->gs_ab[i].dimension();
            const auto dimension_of_ith_ba = this->gs_ba[i].dimension();
            const auto dimension_of_ith_bb = this->gs_bb[i].dimension();
            if ((dimension_of_first_aa != dimension_of_ith_aa) || (dimension_of_first_ab != dimension_of_ith_ab) || (dimension_of_first_ba != dimension_of_ith_ba) || (dimension_of_first_bb != dimension_of_ith_bb)) {
                throw std::invalid_argument("USQTwoElectronOperator(const std::array<QCMatrix<Scalar>, Components>&): The given tenso representations do not have the same dimensions for either the alpha, beta or one of the mixed components.");
            }
        }
    }


    /**
     *  A constructor for ScalarUSQTwoElectronOperators that doesn't require the argument to be a vector of just one element.
     * 
     *  @param g_aa           the tensor representation of the alpha-alpha integrals of this scalar second-quantized operator
     *  @param g_ab           the tensor representation of the alpha-beta integrals of this scalar second-quantized operator
     *  @param g_ba           the tensor representation of the beta-alpha integrals of this scalar second-quantized operator
     *  @param g_bb           the tensor representation of the beta-beta integrals of this scalar second-quantized operator
     * 
     *  @note This constructor is only available for ScalarSQTwoElectronOperators (for the std::enable_if, see https://stackoverflow.com/a/17842695/7930415)
     */
    template <size_t Z = Components>
    USQTwoElectronOperator(const QCRankFourTensor<Scalar>& g_aa, const QCRankFourTensor<Scalar>& g_ab, const QCRankFourTensor<Scalar>& g_ba, const QCRankFourTensor<Scalar>& g_bb, typename std::enable_if<Z == 1>::type* = 0) :
        USQTwoElectronOperator(std::array<QCRankFourTensor<Scalar>, 1>{g_aa}, std::array<QCRankFourTensor<Scalar>, 1>{g_ab}, std::array<QCRankFourTensor<Scalar>, 1>{g_ba}, std::array<QCRankFourTensor<Scalar>, 1>{g_bb})
    {}


    /**
     *  Construct an unrestricted two-electron operator with parameters that are zero. The dimensions of the alpha and beta component are the same.
     * 
     *  @param dim        the dimension of the matrix representation of the alpha and beta parameters, i.e. the number of orbitals/sites
     * 
     */
    USQTwoElectronOperator(const size_t dim) {
        for (size_t i = 0; i < Components; i++) {
            this->gs_aa[i] = QCRankFourTensor<Scalar>(dim);
            this->gs_ab[i] = QCRankFourTensor<Scalar>(dim);
            this->gs_ba[i] = QCRankFourTensor<Scalar>(dim);
            this->gs_bb[i] = QCRankFourTensor<Scalar>(dim);
            this->gs_aa[i].setZero();
            this->gs_ab[i].setZero();
            this->gs_ba[i].setZero();
            this->gs_bb[i].setZero();
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
     *  @param SpinComponent            the requested spin component. The 4 possible combinations can yield all the possible blocks (alpha-alpha, alpha-beta, beta-alpha & beta-beta).
     * 
     *  @return read-only tensor representations of all the parameters (integrals) of the different components of this second-quantized operator, for the requested spin components.
     */
    const std::array<QCRankFourTensor<Scalar>, Components>& allParameters(SpinComponent left, SpinComponent right) const {  

        if (left == SpinComponent::ALPHA && right == SpinComponent::ALPHA) {
            return this->gs_aa;
        }

        else if (left == SpinComponent::ALPHA && right == SpinComponent::BETA) {
            return this->gs_ab;
        }

        else if (left == SpinComponent::BETA && right == SpinComponent::ALPHA) {
            return this->gs_ba;
        }

        else {
            return this->gs_bb;
        };  
    }


    /**
     *  @param SpinComponent            the requested spin component. The 4 possible combinations can yield all the possible blocks (alpha-alpha, alpha-beta, beta-alpha & beta-beta)
     * 
     *  @return the writable tensor representations of all the parameters (integrals) of the different components of this second-quantized operator, for the requested spin components.
     */
    std::array<QCRankFourTensor<Scalar>, Components>& allParameters(SpinComponent left, SpinComponent right) {

        if (left == SpinComponent::ALPHA && right == SpinComponent::ALPHA) {
            return this->gs_aa;
        }

        else if (left == SpinComponent::ALPHA && right == SpinComponent::BETA) {
            return this->gs_ab;
        }

        else if (left == SpinComponent::BETA && right == SpinComponent::ALPHA) {
            return this->gs_ba;
        }

        else {
            return this->gs_bb;
        };
    }

    
    /**
     *  @param d_aa            the alpha-alpha 2-RDM that represents the wave function
     *  @param d_ab            the alpha-beta 2-RDM that represents the wave function
     *  @param d_ba            the beta-alpha 2-RDM that represents the wave function
     *  @param d_bb            the beta-beta 2-RDM that represents the wave function
     *
     *  @return the expectation values of all the components of the two-electron operator, with the given 2-RDMs: this includes the prefactor 1/2
     */
    Vector<Scalar, Components> calculateExpectationValue(const TwoRDM<Scalar>& d_aa, const TwoRDM<Scalar>& d_ab, const TwoRDM<Scalar>& d_ba, const TwoRDM<Scalar>& d_bb) const {

        if ((this->dimension(GQCP::SpinComponent::ALPHA, GQCP::SpinComponent::ALPHA) != d_aa.dimension()) || (this->dimension(GQCP::SpinComponent::ALPHA, GQCP::SpinComponent::ALPHA) != d_ab.dimension()) || (this->dimension(GQCP::SpinComponent::ALPHA, GQCP::SpinComponent::ALPHA) != d_ba.dimension()) || (this->dimension(GQCP::SpinComponent::ALPHA, GQCP::SpinComponent::ALPHA) != d_bb.dimension())) {
            throw std::invalid_argument("USQTwoElectronOperator::calculateExpectationValue(const TwoRDM<double>&): One of the given 2-RDMs is not compatible with the respective component of the two-electron operator.");
        }


        std::array<Scalar, Components> expectation_values {};
        for (size_t i = 0; i < Components; i++) {

            // Specify the contractions for the relevant contraction of the two-electron integrals and the 2-RDMs
            //      0.5 g(p q r s) d(p q r s)
            Eigen::array<Eigen::IndexPair<int>, 4> contractions = {Eigen::IndexPair<int>(0,0), Eigen::IndexPair<int>(1,1), Eigen::IndexPair<int>(2,2), Eigen::IndexPair<int>(3,3)};
            //      Perform the contraction
            Eigen::Tensor<Scalar, 0> contraction_aa = 0.5 * this->parameters(GQCP::SpinComponent::ALPHA, GQCP::SpinComponent::ALPHA, i).contract(d_aa.Eigen(), contractions);
            Eigen::Tensor<Scalar, 0> contraction_ab = 0.5 * this->parameters(GQCP::SpinComponent::ALPHA, GQCP::SpinComponent::BETA, i).contract(d_ab.Eigen(), contractions);
            Eigen::Tensor<Scalar, 0> contraction_ba = 0.5 * this->parameters(GQCP::SpinComponent::BETA, GQCP::SpinComponent::ALPHA, i).contract(d_ba.Eigen(), contractions);
            Eigen::Tensor<Scalar, 0> contraction_bb = 0.5 * this->parameters(GQCP::SpinComponent::BETA, GQCP::SpinComponent::BETA, i).contract(d_bb.Eigen(), contractions);

            // As the contraction is a scalar (a tensor of rank 0), we should access by (0).
            expectation_values[i] = contraction_aa(0) + contraction_ab(0) + contraction_ba(0) + contraction_bb(0);
        }

        return Eigen::Map<Eigen::Matrix<Scalar, Components, 1>>(expectation_values.data());  // convert std::array to Vector
    }


    /**
     *  @param SpinComponent            the requested spin component. The 4 possible combinations can yield all the possible blocks (alpha-alpha, alpha-beta, beta-alpha & beta-beta)
     * 
     *  @return the dimension of the tensors for the requested spin components.
     */
    size_t dimension(SpinComponent left, SpinComponent right) const {

         if (left == SpinComponent::ALPHA && right == SpinComponent::ALPHA) {
            return this->gs_aa[0].dimension();
        }

        else if (left == SpinComponent::ALPHA && right == SpinComponent::BETA) {
            return this->gs_ab[0].dimension();
        }

        else if (left == SpinComponent::BETA && right == SpinComponent::ALPHA) {
            return this->gs_ba[0].dimension();
        }

        else {
            return this->gs_bb[0].dimension();
        };
    }


    /**
     *  @param i                        The index of the component.
     *  @param SpinComponent            the requested spin component. The 4 possible combinations can yield all the possible blocks (alpha-alpha, alpha-beta, beta-alpha & beta-beta)
     * 
     *  @return a read-only tensor representation of the parameters (integrals) of one of the the different components of this second-quantized operator, for the requested spin components.
     */
    const QCRankFourTensor<Scalar>& parameters(SpinComponent left, SpinComponent right, const size_t i = 0 ) const {

        if (left == SpinComponent::ALPHA && right == SpinComponent::ALPHA) {
            return this->gs_aa[i];
        }

        else if (left == SpinComponent::ALPHA && right == SpinComponent::BETA) {
            return this->gs_ab[i];
        }

        else if (left == SpinComponent::BETA && right == SpinComponent::ALPHA) {
            return this->gs_ba[i];
        }

        else {
            return this->gs_bb[i];
        };
    }


    /**
     *  @param i                        The index of the component.
     *  @param SpinComponent            the requested spin component. The 4 possible combinations can yield all the possible blocks (alpha-alpha, alpha-beta, beta-alpha & beta-beta)
     * 
     *  @return a writable tensor representation of the parameters (integrals) of one of the the different components of this second-quantized operator, for the requested spin components.
     */
    QCRankFourTensor<Scalar>& parameters(SpinComponent left, SpinComponent right, const size_t i = 0) {

        if (left == SpinComponent::ALPHA && right == SpinComponent::ALPHA) {
            return this->gs_aa[i];
        }

        else if (left == SpinComponent::ALPHA && right == SpinComponent::BETA) {
            return this->gs_ab[i];
        }

        else if (left == SpinComponent::BETA && right == SpinComponent::ALPHA) {
            return this->gs_ba[i];
        }

        else {
            return this->gs_bb[i];
        };
    }


    /**
     *  In-place rotate the operator to another basis. The alpha-alpha, alpha-beta, beta-alpha and beta-beta components are transformed in the same way.
     * 
     *  @param U                            the (unitary) rotation matrix
     */
    void rotate(const TransformationMatrix<Scalar>& U) {

        // Transform the matrix representations of the components
        for (auto& g_aa : this->allParameters(GQCP::SpinComponent::ALPHA, GQCP::SpinComponent::ALPHA)) {
            g_aa.basisRotateInPlace(U);
        }
        for (auto& g_ab : this->allParameters(GQCP::SpinComponent::ALPHA, GQCP::SpinComponent::BETA)) {
            g_ab.basisRotateInPlace(U);
        }
        for (auto& g_ba : this->allParameters(GQCP::SpinComponent::BETA, GQCP::SpinComponent::ALPHA)) {
            g_ba.basisRotateInPlace(U);
        }
        for (auto& g_bb : this->allParameters(GQCP::SpinComponent::BETA, GQCP::SpinComponent::BETA)) {
            g_bb.basisRotateInPlace(U);
        }
    }


    /**
     *  In-place rotate the operator using a unitary Jacobi rotation matrix constructed from the Jacobi rotation parameters. The alpha-alpha, alpha-beta, beta-alpha and beta-beta components are transformed in the same way.
     * 
     *  @param jacobi_rotation_parameters       the Jacobi rotation parameters (p, q, angle) that are used to specify a Jacobi rotation: we use the (cos, sin, -sin, cos) definition for the Jacobi rotation matrix
     */
    void rotate(const JacobiRotationParameters& jacobi_rotation_parameters) {

        // Transform the matrix representations of the components
        for (auto& g_aa : this->allParameters(GQCP::SpinComponent::ALPHA, GQCP::SpinComponent::ALPHA)) {
            g_aa.basisRotateInPlace(jacobi_rotation_parameters);
        }
        for (auto& g_ab : this->allParameters(GQCP::SpinComponent::ALPHA, GQCP::SpinComponent::BETA)) {
            g_ab.basisRotateInPlace(jacobi_rotation_parameters);
        }
        for (auto& g_ba : this->allParameters(GQCP::SpinComponent::BETA, GQCP::SpinComponent::ALPHA)) {
            g_ba.basisRotateInPlace(jacobi_rotation_parameters);
        }
        for (auto& g_bb : this->allParameters(GQCP::SpinComponent::BETA, GQCP::SpinComponent::BETA)) {
            g_bb.basisRotateInPlace(jacobi_rotation_parameters);
        }
    }


    /**
     *  In-place transform the operator to another basis. The alpha-alpha, alpha-beta, beta-alpha and beta-beta components are transformed in the same way.
     * 
     *  @param T                            the transformation matrix
     */
    void transform(const TransformationMatrix<Scalar>& T) {

        // Transform the matrix representations of the components
        for (auto& g_aa : this->allParameters(GQCP::SpinComponent::ALPHA, GQCP::SpinComponent::ALPHA)) {
            g_aa.basisTransformInPlace(T);
        }
        for (auto& g_ab : this->allParameters(GQCP::SpinComponent::ALPHA, GQCP::SpinComponent::BETA)) {
            g_ab.basisTransformInPlace(T);
        }
        for (auto& g_ba : this->allParameters(GQCP::SpinComponent::BETA, GQCP::SpinComponent::ALPHA)) {
            g_ba.basisTransformInPlace(T);
        }
        for (auto& g_bb : this->allParameters(GQCP::SpinComponent::BETA, GQCP::SpinComponent::BETA)) {
            g_bb.basisTransformInPlace(T);
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

    auto G_sum_aa = lhs.allParameters(GQCP::SpinComponent::ALPHA, GQCP::SpinComponent::ALPHA);
    auto G_sum_ab = lhs.allParameters(GQCP::SpinComponent::ALPHA, GQCP::SpinComponent::BETA);
    auto G_sum_ba = lhs.allParameters(GQCP::SpinComponent::BETA, GQCP::SpinComponent::ALPHA);
    auto G_sum_bb = lhs.allParameters(GQCP::SpinComponent::BETA, GQCP::SpinComponent::BETA);
    for (size_t i = 0; i < Components; i++) {
        G_sum_aa[i] += rhs.parameters(GQCP::SpinComponent::ALPHA, GQCP::SpinComponent::ALPHA, i);
        G_sum_ab[i] += rhs.parameters(GQCP::SpinComponent::ALPHA, GQCP::SpinComponent::BETA, i);
        G_sum_ba[i] += rhs.parameters(GQCP::SpinComponent::BETA, GQCP::SpinComponent::ALPHA, i);
        G_sum_bb[i] += rhs.parameters(GQCP::SpinComponent::BETA, GQCP::SpinComponent::BETA, i);
    }

    return USQTwoElectronOperator<ResultScalar, Components>(G_sum_aa, G_sum_ab, G_sum_ba, G_sum_bb);
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

    auto G_aa = op.allParameters(GQCP::SpinComponent::ALPHA, GQCP::SpinComponent::ALPHA);
    auto G_ab = op.allParameters(GQCP::SpinComponent::ALPHA, GQCP::SpinComponent::BETA);
    auto G_ba = op.allParameters(GQCP::SpinComponent::BETA, GQCP::SpinComponent::ALPHA);
    auto G_bb = op.allParameters(GQCP::SpinComponent::BETA, GQCP::SpinComponent::BETA);
    for (size_t i = 0; i < Components; i++) {
        G_aa[i] = scalar * G_aa[i];
        G_ab[i] = scalar * G_ab[i];
        G_ba[i] = scalar * G_ba[i];
        G_bb[i] = scalar * G_bb[i];
    }

    return USQTwoElectronOperator<ResultScalar, Components>(G_aa, G_ab, G_ba, G_bb);
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
