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
 *  A class that represents an unrestriccted second-quantized two-electron operator: it holds the alpha and beta matrix representation of its parameters, which are (usually) integrals over first-quantized operators
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
    std::array<QCRankFourTensor<Scalar>, Components> gs_alpha;  // all the matrix representations (hence the 's') of the alpha parameters (integrals) of the different components of this second-quantized operator
    std::array<QCRankFourTensor<Scalar>, Components> gs_beta;  // all the matrix representations (hence the 's') of the beta parameters (integrals) of the different components of this second-quantized operator

public:

    /*
     *  CONSTRUCTORS
     */
    /**
     *  @param gs_alpha            all the matrix representations (hence the 's') of the alpha parameters (integrals) of the different components of this second-quantized operator
     *  @param gs_beta             all the matrix representations (hence the 's') of the beta parameters (integrals) of the different components of this second-quantized operator
     * 
     */
    USQTwoElectronOperator(const std::array<QCRankFourTensor<Scalar>, Components>& gs_alpha, const std::array<QCRankFourTensor<Scalar>, Components>& gs_beta) : 
        gs_alpha (gs_alpha),
        gs_beta (gs_beta)
    {
        // Check if the given matrix representations have the same dimensions
        const auto dimension_of_first_alpha = this->gs_alpha[0].dimension();
        const auto dimension_of_first_beta = this-> gs_beta[0].dimension();

        for (size_t i = 1; i < Components; i++) {

            const auto dimension_of_ith_alpha = this->gs_alpha[i].dimension();
            const auto dimension_of_ith_beta = this->gs_beta[i].dimension();
            if ((dimension_of_first_alpha != dimension_of_ith_alpha) || (dimension_of_first_beta != dimension_of_ith_beta)) {
                throw std::invalid_argument("SQTwoElectronOperator(const std::array<QCMatrix<Scalar>, Components>&): The given matrix representations do not have the same dimensions for either the alpha or beta component.");
            }
        }
    }


    /**
     *  A constructor for ScalarUSQTwoElectronOperators that doesn't require the argument to be a vector of just one element.
     * 
     *  @param g_alpha           the matrix representation of the alpha integrals of this scalar second-quantized operator
     *  @param g_beta            the matrix representation of the beta integrals of this scalar second-quantized operator
     * 
     *  @note This constructor is only available for ScalarSQTwoElectronOperators (for the std::enable_if, see https://stackoverflow.com/a/17842695/7930415)
     */
    template <size_t Z = Components>
    USQTwoElectronOperator(const QCRankFourTensor<Scalar>& g_alpha, const QCRankFourTensor<Scalar>& g_beta, typename std::enable_if<Z == 1>::type* = 0) :
        USQTwoElectronOperator(std::array<QCRankFourTensor<Scalar>, 1>{g_alpha}, std::array<QCRankFourTensor<Scalar>, 1>{g_beta})
    {}


    /**
     *  Construct an unrestricted two-electron operator with zero parameters, dimensions of alpha and beta component may vary.
     * 
     *  @param dim_alpha          the dimension of the matrix representation of the alpha parameters, i.e. the number of orbitals/sites
     *  @param dim_beta           the dimension of the matrix representation of the beta parameters, i.e. the number of orbitals/sites

     */
    USQTwoElectronOperator(const size_t dim_alpha, const size_t dim_beta) {
        for (size_t i = 0; i < Components; i++) {
            this->gs_alpha[i] = QCRankFourTensor<Scalar>(dim_alpha);
            this->gs_beta[i] = QCRankFourTensor<Scalar>(dim_beta);
            this->gs_alpha[i].setZero();
            this->gs_beta[i].setZero();
        }
    }


    /**
     *  Construct an unrestricted two-electron operator with zero parameters, dimensions of alpha and beta component are the same.
     * 
     *  @param dim        the dimension of the matrix representation of the alpha and beta parameters, i.e. the number of orbitals/sites
     * 
     */
    USQTwoElectronOperator(const size_t dim) :
        USQTwoElectronOperator(dim, dim)
        {}
    

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
     *  @return read-only matrix representations of all the alpha parameters (integrals) of the different components of this second-quantized operator
     */
    const std::array<QCRankFourTensor<Scalar>, Components>& allAlphaParameters() const {
        return this->gs_alpha;
    }

    /**
     *  @return writable matrix representations of all the alpha parameters (integrals) of the different components of this second-quantized operator
     */
    std::array<QCRankFourTensor<Scalar>, Components>& allAlphaParameters() {
        return this->gs_alpha;
    }


    /**
     *  @return read-only matrix representations of all the beta parameters (integrals) of the different components of this second-quantized operator
     */
    const std::array<QCRankFourTensor<Scalar>, Components>& allBetaParameters() const {
        return this->gs_beta;
    }


    /**
     *  @return writable matrix representations of all the beta parameters (integrals) of the different components of this second-quantized operator
     */
    std::array<QCRankFourTensor<Scalar>, Components>& allBetaParameters() {
        return this->gs_beta;
    }


    /**
     *  @return the dimension of the alpha components
     */
    size_t alphaDimension() const { return this->gs_alpha[0].dimension(); }


    /**
     *  @param i            the index of the component
     * 
     *  @return a read-only matrix representation of the alpha parameters (integrals) of one of the the different components of this second-quantized operator
     */
    const QCRankFourTensor<Scalar>& alphaParameters(const size_t i = 0) const {
        return this->gs_alpha[i];
    }


    /**
     *  @param i            the index of the component
     * 
     *  @return a writable the matrix representation of the alpha parameters (integrals) of one of the the different components of this second-quantized operator
     */
    QCRankFourTensor<Scalar>& alphaParameters(const size_t i = 0) {
        return this->gs_alpha[i];
    }


    /**
     *  @return the dimension of the alpha components
     */
    size_t betaDimension() const { return this->gs_beta[0].dimension(); }


    /**
     *  @param i            the index of the component
     * 
     *  @return a read-only matrix representation of the beta parameters (integrals) of one of the the different components of this second-quantized operator
     */
    const QCRankFourTensor<Scalar>& betaParameters(const size_t i = 0) const {
        return this->gs_beta[i];
    }


    /**
     *  @param i            the index of the component
     * 
     *  @return a writable the matrix representation of the beta parameters (integrals) of one of the the different components of this second-quantized operator
     */
    QCRankFourTensor<Scalar>& betaParameters(const size_t i = 0) {
        return this->gs_beta[i];
    }

    
    /**
     *  @param d            the 2-RDM that represents the wave function
     *
     *  @return the expectation values of all the components of the two-electron operator, with the given 2-RDM: this includes the prefactor 1/2
     */
    Vector<Scalar, Components> calculateExpectationValue(const TwoRDM<Scalar>& d_alpha, const TwoRDM<Scalar>& d_beta) const {

        if ((this->alphaDimension() != d_alpha.dimension()) || this->betaDimension() != d_beta.dimension()) {
            throw std::invalid_argument("SQTwoElectronOperator::calculateExpectationValue(const TwoRDM<double>&): The given alpha or beta 2-RDM is not compatible with the respective two-electron operator.");
        }


        std::array<Scalar, Components> expectation_values {};  // zero initialization of alphha component
        for (size_t i = 0; i < Components; i++) {

            // Specify the contractions for the relevant contraction of the two-electron integrals and the 2-RDM
            //      0.5 g(p q r s) d(p q r s)
            Eigen::array<Eigen::IndexPair<int>, 4> contractions = {Eigen::IndexPair<int>(0,0), Eigen::IndexPair<int>(1,1), Eigen::IndexPair<int>(2,2), Eigen::IndexPair<int>(3,3)};
            //      Perform the contraction
            Eigen::Tensor<Scalar, 0> contraction_alpha = 0.5 * this->alphaParameters(i).contract(d_alpha.Eigen(), contractions);
            Eigen::Tensor<Scalar, 0> contractions_beta = 0.5 * this->betaParameters(i).contract(d_beta.Eigen(), contractions);

            // As the contraction is a scalar (a tensor of rank 0), we should access by (0).
            expectation_values[i] = contraction_alpha(0) + contractions_beta(0);
        }

        return Eigen::Map<Eigen::Matrix<Scalar, Components, 1>>(expectation_values.data());  // convert std::array to Vector
    }


    /**
     *  @return the one-electron operator that is the difference between a two-electron operator (e_pqrs) and a product of one-electron operators (E_pq E_rs)
     */
    USQOneElectronOperator<Scalar, Components> effectiveOneElectronPartition() const {

        // Initialize a zero operator
        const auto K_alpha = this->alphaDimension();  // number of alpha orbitals
        const auto K_beta = this->betaDimension(); // number of beta orbitals
        USQOneElectronOperator<Scalar, Components> F (K_alpha, K_beta);


        // Use a formula to set the parameters
        for (size_t i = 0; i < Components; i++) {

            for (size_t p = 0; p < K_alpha; p++) {
                for (size_t q = 0; q < K_alpha; q++) {
                    for (size_t r = 0; r < K_alpha; r++) {
                        F.alphaParameters(i)(p,q) -= 0.5 * this->alphaParameters(i)(p,r,r,q);
                        F.betaParameters(i)(p,q) -= 0.5 * this->betaParameters(i)(p,r,r,q);
                    }
                }
            }

        }  // loop over components

        return F;
    }


    /**
     *  In-place rotate the operator to another basis. The alpha and beta components are transformed in the same way.
     * 
     *  @param U                            the (unitary) rotation matrix
     */
    void rotate(const TransformationMatrix<Scalar>& U) {

        // Transform the matrix representations of the components
        for (auto& g_alpha : this->allAlphaParameters()) {
            g_alpha.basisRotateInPlace(U);
        }
        for (auto& g_beta : this->allBetaParameters()) {
            g_beta.basisRotateInPlace(U);
        }
    }


    /**
     *  In-place rotate the operator using a unitary Jacobi rotation matrix constructed from the Jacobi rotation parameters. The alpha and beta components are transformed in the same way.
     * 
     *  @param jacobi_rotation_parameters       the Jacobi rotation parameters (p, q, angle) that are used to specify a Jacobi rotation: we use the (cos, sin, -sin, cos) definition for the Jacobi rotation matrix
     */
    void rotate(const JacobiRotationParameters& jacobi_rotation_parameters) {

        // Transform the matrix representations of the components
        for (auto& g_alpha : this->allAlphaParameters()) {
            g_alpha.basisRotateInPlace(jacobi_rotation_parameters);
        }
        for (auto& g_beta : this->allBetaParameters()) {
            g_beta.basisRotateInPlace(jacobi_rotation_parameters);
        }
    }


    /**
     *  In-place transform the operator to another basis. The alpha and beta components are transformed in the same way.
     * 
     *  @param T                            the transformation matrix
     */
    void transform(const TransformationMatrix<Scalar>& T) {

        // Transform the matrix representations of the components
        for (auto& g_alpha : this->allAlphaParameters()) {
            g_alpha.basisTransformInPlace(T);
        }
        for (auto& g_beta : this->allBetaParameters()) {
            g_beta.basisTransformInPlace(T);
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
 *  Add two two-electron operators by adding their parameters. The two alpha components are added together and the two beta components are added together.
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

    auto G_sum_alpha = lhs.allAlphaParameters();
    auto G_sum_beta = lhs.allBetaParameters();
    for (size_t i = 0; i < Components; i++) {
        G_sum_alpha[i] += rhs.alphaParameters(i);
        G_sum_beta[i] += rhs.betaParameters(i);
    }

    return USQTwoElectronOperator<ResultScalar, Components>(G_sum_alpha, G_sum_beta);
}


/**
 *  Multiply a two-electron operator with a scalar. The alpha and beta component are multiplied with the same scalar.
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

    auto G_alpha = op.allAlphaParameters();
    auto G_beta = op.allBetaParameters();
    for (size_t i = 0; i < Components; i++) {
        G_alpha[i] = scalar * G_alpha[i];
        G_beta[i] = scalar * G_beta[i];
    }

    return USQTwoElectronOperator<ResultScalar, Components>(G_alpha, G_beta);
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
 *  Subtract two two-electron operators by subtracting their respective alpha and beta parameters
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
