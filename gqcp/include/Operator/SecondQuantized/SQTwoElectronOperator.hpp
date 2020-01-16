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
#include "Operator/SecondQuantized/SQOneElectronOperator.hpp"
#include "Processing/RDM/OneRDM.hpp"
#include "Processing/RDM/TwoRDM.hpp"
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
    std::array<QCRankFourTensor<Scalar>, Components> gs;  // all the matrix representations (hence the 's') of the parameters (integrals) of the different components of this second-quantized operator

public:

    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param gs            all the matrix representations (hence the 's') of the parameters (integrals) of the different components of this second-quantized operator
     */
    SQTwoElectronOperator(const std::array<QCRankFourTensor<Scalar>, Components>& gs) : 
        gs (gs)
    {
        // Check if the given matrix representations have the same dimensions
        const auto dimension_of_first = this->gs[0].dimension();

        for (size_t i = 1; i < Components; i++) {

            const auto dimension_of_ith = this->gs[i].dimension();
            if (dimension_of_first != dimension_of_ith) {
                throw std::invalid_argument("SQTwoElectronOperator(const std::array<QCMatrix<Scalar>, Components>&): The given matrix representations do not have the same dimensions.");
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
            this->gs[i] = QCRankFourTensor<Scalar>(dim);
            this->gs[i].setZero();
        }
    }


    /**
     *  Default constructor: construct a two-electron operator with parameters that are zero
     */
    SQTwoElectronOperator() :
        SQTwoElectronOperator(0)  // dimensions of the representations are zero
    {}



    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return the dimension of the matrix representation of the parameters, i.e. the number of orbitals/sites
     */
    size_t dimension() const {
        return this->gs[0].dimension();  // all dimensions are the same, this is checked in the constructors
    }

    size_t get_dim() const {
        return this->dimension();
    }

    size_t get_K() const {
        return this->dimension();
    }


    /**
     *  @return read-only matrix representations of all the parameters (integrals) of the different components of this second-quantized operator
     */
    const std::array<QCRankFourTensor<Scalar>, Components>& allParameters() const {
        return this->gs;
    }


    /**
     *  @return writable matrix representations of all the parameters (integrals) of the different components of this second-quantized operator
     */
    std::array<QCRankFourTensor<Scalar>, Components>& allParameters() {
        return this->gs;
    }


    /**
     *  @param i            the index of the component
     * 
     *  @return a read-only matrix representation of the parameters (integrals) of one of the the different components of this second-quantized operator
     */
    const QCRankFourTensor<Scalar>& parameters(const size_t i = 0) const {
        return this->gs[i];
    }


    /**
     *  @param i            the index of the component
     * 
     *  @return a writable the matrix representation of the parameters (integrals) of one of the the different components of this second-quantized operator
     */
    QCRankFourTensor<Scalar>& parameters(const size_t i = 0) {
        return this->gs[i];
    }


    /**
     *  In-place transform the operator to another basis
     * 
     *  @param T                            the transformation matrix
     */
    void transform(const TransformationMatrix<Scalar>& T) {

        // Transform the matrix representations of the components
        for (auto& g : this->allParameters()) {
            g.basisTransformInPlace(T);
        }
    }


    /**
     *  In-place rotate the operator to another basis
     * 
     *  @param U                            the (unitary) rotation matrix
     */
    void rotate(const TransformationMatrix<Scalar>& U) {

        // Transform the matrix representations of the components
        for (auto& g : this->allParameters()) {
            g.basisRotateInPlace(U);
        }
    }


    /**
     *  In-place rotate the operator using a unitary Jacobi rotation matrix constructed from the Jacobi rotation parameters
     * 
     *  @param jacobi_rotation_parameters       the Jacobi rotation parameters (p, q, angle) that are used to specify a Jacobi rotation: we use the (cos, sin, -sin, cos) definition for the Jacobi rotation matrix
     */
    void rotate(const JacobiRotationParameters& jacobi_rotation_parameters) {

        // Transform the matrix representations of the components
        for (auto& g : this->allParameters()) {
            g.basisRotateInPlace(jacobi_rotation_parameters);
        }
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
                        F.parameters(i)(p,q) -= 0.5 * this->parameters(i)(p,r,r,q);
                    }
                }
            }

        }  // loop over components

        return F;
    }


    /**
     *  @param D      the 1-DM (or the response 1-DM for made-variational wave function models)
     *  @param d      the 2-DM (or the response 2-DM for made-variational wave function models)
     *
     *  @return the (generalized) Fockian matrix for each of the components
     */
    std::array<SquareMatrix<Scalar>, Components> calculateFockianMatrix(const OneRDM<double>& D, const TwoRDM<double>& d) const {

        // Check if dimensions are compatible
        if (D.dimension() != this->dimension()) {
            throw std::invalid_argument("SQTwoElectronOperator::calculateFockianMatrix(OneRDM<double>, TwoRDM<double>): The 1-RDM is not compatible with the two-electron operator.");
        }

        if (d.dimension() != this->dimension()) {
            throw std::invalid_argument("SQTwoElectronOperator::calculateFockianMatrix(OneRDM<double>, TwoRDM<double>): The 2-RDM is not compatible with the two-electron operator.");
        }


        // A KISS implementation of the calculation of the generalized Fock matrix F
        std::array<SquareMatrix<Scalar>, Components> Fs;  // Fock matrices (hence the 's')
        for (size_t i = 0; i < Components; i++) {

            const auto& g_i = this->parameters(i);  // the matrix representation of the parameters of the i-th component

            // Calculate the Fockian matrix for every component and add it to the array
            SquareMatrix<Scalar> F_i = SquareMatrix<Scalar>::Zero(this->dimension(), this->dimension());
            for (size_t p = 0; p < this->dimension(); p++) {
                for (size_t q = 0; q < this->dimension(); q++) {

                    for (size_t r = 0; r < this->dimension(); r++) {
                        for (size_t s = 0; s < this->dimension(); s++) {
                            for (size_t t = 0; t < this->dimension(); t++) {
                                F_i(p,q) += g_i(q,r,s,t) * (d(p,r,s,t) + d(r,p,s,t));
                            }
                        }
                    }

                }
            }  // F elements loop
            Fs[i] = 0.5 * F_i;
        }

        return Fs;
    }


    /**
     *  @param D      the 1-DM (or the response 1-DM for made-variational wave function models)
     *  @param d      the 2-DM (or the response 2-DM for made-variational wave function models)
     *
     *  @return the (generalized) super-Fockian matrix
     */
    std::array<SquareRankFourTensor<Scalar>, Components> calculateSuperFockianMatrix(const OneRDM<double>& D, const TwoRDM<double>& d) const {

        // Check if dimensions are compatible
        if (D.dimension() != this->dimension()) {
            throw std::invalid_argument("SQOneElectronOperator::calculateFockianMatrix(OneRDM<double>, TwoRDM<double>): The 1-RDM is not compatible with the one-electron operator.");
        }

        if (d.dimension() != this->dimension()) {
            throw std::invalid_argument("SQOneElectronOperator::calculateFockianMatrix(OneRDM<double>, TwoRDM<double>): The 2-RDM is not compatible with the one-electron operator.");
        }


        // A KISS implementation of the calculation of the super-Fockian matrix
        std::array<SquareRankFourTensor<Scalar>, Components> Gs;  // multiple Gs, hence the 's'
        const auto Fs = this->calculateFockianMatrix(D, d);  // the Fockian matrices are necessary in the calculation
        for (size_t i = 0; i < Components; i++) {
            
            const auto& g_i = this->parameters(i);  // the matrix representation of the parameters of the i-th component
            const auto& F_i = Fs[i];  // the Fockian matrix of the i-th component

            // Calculate the super-Fockian matrix for every component and add it to the array
            SquareRankFourTensor<Scalar> G_i (this->dimension());
            G_i.setZero();
            for (size_t p = 0; p < this->dimension(); p++) {
                for (size_t q = 0; q < this->dimension(); q++) {
                    for (size_t r = 0; r < this->dimension(); r++) {
                        for (size_t s = 0; s < this->dimension(); s++) {

                            if (q == r) {
                                G_i(p,q,r,s) += 2 * F_i(p,s);
                            }

                            for (size_t t = 0; t < this->dimension(); t++) {
                                for (size_t u = 0; u < this->dimension(); u++) {
                                    G_i(p,q,r,s) += g_i(s,t,q,u) * (d(r,t,p,u) + d(t,r,u,p)) - g_i(s,t,u,p) * (d(r,t,u,q) + d(t,r,q,u)) - g_i(s,p,t,u) * (d(r,q,t,u) + d(q,r,u,t));
                                }
                            }

                        }
                    }
                }
            }  // G_i elements loop
            Gs[i] = 0.5 * G_i;
        }

        return Gs;
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


/**
 *  Negate a two-electron operator
 * 
 *  @tparam Scalar              the scalar type of the operator
 *  @tparam Components          the number of components of the one-electron operator
 * 
 *  @param op                   the operator
 */
template <typename Scalar, size_t Components>
SQTwoElectronOperator<Scalar, Components> operator-(const SQTwoElectronOperator<Scalar, Components>& op) {

    return (-1.0) * op;  // negation is scalar multiplication with (-1.0)
}


/**
 *  Subtract two two-electron operators by adding their parameters
 * 
 *  @tparam LHSScalar           the scalar type of the left-hand side
 *  @tparam RHSScalar           the scalar type of the right-hand side
 *  @tparam Components          the number of components of the two-electron operators
 * 
 *  @param lhs                  the left-hand side
 *  @param rhs                  the right-hand side
 */
template <typename LHSScalar, typename RHSScalar, size_t Components>
auto operator-(const SQTwoElectronOperator<LHSScalar, Components>& lhs, const SQTwoElectronOperator<RHSScalar, Components>& rhs) -> SQTwoElectronOperator<sum_t<LHSScalar, RHSScalar>, Components> {

    return lhs + (-rhs);
}


}  // namespace GQCP
