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
#include "Mathematical/Representation/QCMatrix.hpp"
#include "Mathematical/ScalarFunction.hpp"
#include "Processing/RDM/OneRDM.hpp"
#include "Processing/RDM/TwoRDM.hpp"
#include "Utilities/type_traits.hpp"
#include "Basis/SpinorBasis/SpinComponent.hpp"

#include <array>


namespace GQCP {


/**
 *  A class that represents an unrestricted second-quantized one-electron operator: it holds the matrix representation of its parameters for both spin components, which are (usually) integrals over first-quantized operators
 *
 *  @tparam _Scalar             the scalar type, i.e. the scalar representation of one of the parameters
 *  @tparam _Components         the number of components of the second-quantized operator
 */
template <typename _Scalar, size_t _Components, _Scalar, size_t _Components>
class USQOneElectronOperator{
public:

    using Scalar = _Scalar;
    static constexpr auto Components = _Components;

private:
    std::array<QCMatrix<Scalar>, Components> fs_alpha;  // all the matrix representations (hence the s) of the parameters (integrals) of the different components of this second-quantized operator
    std::array<QCMatrix<Scalar>, Components> fs_beta;  // all the matrix representations (hence the s) of the parameters (integrals) of the different components of this second-quantized operator

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
        // Check if the given matrix representations have the same dimensions
        const auto dimension_of_first_alpha = this->fs_alpha[0].dimension();
        const auto dimension_of_first_beta = this->fs_beta[0].dimension();
        for (size_t i = 1; i < Components; i++) {

            const auto dimension_of_ith_alpha = this->fs_alpha[i].dimension();
            const auto dimension_of_ith_beta = this->fs_beta[i].dimension();

            if (dimension_of_first_alpha != dimension_of_ith_alpha || dimension_of_first_beta != dimension_of_ith_beta){
                throw std::invalid_argument("USQOneElectronOperator(const std::array<QCMatrix<Scalar>, Components>&, const std::array<QCMatrix<Scalar>, components>&): 
                The given matrix representations do not have the same dimensions for either the alpha or beta component.");
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
    USQOneElectronOperator(const QCMatrix<Scalar>& f_alpha, typename std::enable_if<Z == 1>::type* = 0, const QCMatrix<Scalar>& f_beta, typename std::enable_if<Z == 1>::type* = 0) :
        USQOneElectronOperator(std::array<QCMatrix<Scalar>, 1>{f_alpha}, std::array<QCMatrix<Scalar>, 1>{f_beta})
    {}


    /**
     *  Construct a one-electron operator with parameters that are zero, for both the alpha and beta spin components
     * 
     *  @param dim          the dimension of the matrix representation of the parameters, i.e. the number of orbitals/sites. This dimension is the same for the alpha and beta component.
     */
    USQOneElectronOperator(const size_t dim, const size_t dim) {
        for (size_t i = 0; i < Components; i++) {
            this->fs_alpha[i] = QCMatrix<Scalar>::Zero(dim, dim);
            this->fs_beta[i] = QCMatrix<Scalar>::Zero(dim, dim);
        }
    }


    /**
     *  Default constructor: construct a one-electron operator with parameters that are zero, for both spin components.
     */
    USQOneElectronOperator() :
        USQOneElectronOperator(0, 0)  // dimensions of the representations are zero
    {}


    /*
     *  GETTERS
     */
    size_t get_dim_alpha() const { return this->fs_alpha.dimension(); }
    size_t get_dim_beta() const { return this->fs_beta.dimension(); }
    size_t get_K_alpha() const { return this->fs_alpha.dimension(); }
    size_t get_K_beta() const { return this->fs_beta.dimension(); }


    /*
     *  OPERATORS
     */
    /**
     *  @param i            the index
     * 
     *  @return the i-th component of this operator, for both spin components
     */
    USQOneElectronOperator<Scalar, 1, Scalar, 1> operator[](const size_t i) const {

        if (i >= Components) {
            throw std::invalid_argument("USQOneElectronOperator::operator[](const size_t): The given index is out of bounds.");
        }

        return USQOneElectronOperator<Scalar, 1, Scalar, 1> {this->fs_alpha[i], this->fs_beta[i]};
    }


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return read-only matrix representations of all the parameters (integrals) of the different components of this second-quantized operator
     */
    const std::array<QCMatrix<Scalar>, Components>& allParameters() const {
        return this->fs_alpha, this->fs_beta;
    }


    /**
     *  @return writable matrix representations of all the parameters (integrals) of the different components of this second-quantized operator
     */
    std::array<QCMatrix<Scalar>, Components>& allParameters() {
        return this->fs_alpha, this->fs_beta;
    }


    /**
     *  @param D                the 1-RDM that represents the wave function
     *
     *  @return the expectation values of all components of the one-electron operator
     */
    Vector<Scalar, Components> calculateExpectationValue(const OneRDM<Scalar>& D_alpha, const OneRDM<Scalar>& D_beta) const {

        if (this->fs_alpha.dimension() != D_alpha.dimension() || this->fs_beta.dimension() != D_beta.dimension()) {
            throw std::invalid_argument("USQOneElectronOperator::calculateExpectationValue(const OneRDM<Scalar>, OneRDM<Scalar>): The given 1-RDM is not compatible with the one-electron operator.");
        }

        std::array<Scalar, Components> expectation_values {};  // zero initialization
        for (size_t i = 0; i < Components; i++) {
            expectation_values_alpha[i] = (this->fs_alpha.parameters(i) * D_alpha).trace();
            expectation_values_beta[i] = (this->fs_beta.parameters(i) * D_beta).trace();
        }

        return Eigen::Map<Eigen::Matrix<Scalar, Components, 1>>(expectation_values_alpha.data()), Eigen::Map<Eigen::Matrix<Scalar, Components, 1>>(expectation_values_beta.data());  // convert std::array to Vector
    }


    /**
     *  @param D      the 1-DM (or the response 1-DM for made-variational wave function models)
     *  @param d      the 2-DM (or the response 2-DM for made-variational wave function models)
     *
     *  @return the (generalized) Fockian matrix for each of the components
     */
    std::array<SquareMatrix<Scalar>, Components> calculateFockianMatrix(const OneRDM<double>& D_alpha, const TwoRDM<double>& d_alpha, const OneRDM<double>& D_beta, const TwoRDM<double>& d_beta) const {

        // Check if dimensions are compatible
        if (D_alpha.dimension() != this->fs_alpha.dimension() || D_beta.dimension() != this->fs_beta.dimension()) {
            throw std::invalid_argument("USQOneElectronOperator::calculateFockianMatrix(OneRDM<double>, TwoRDM<double>, OneRDM<double>, TwoRDM<double>): The 1-RDM is not compatible with the one-electron operator.");
        }

        if (d_alpha.dimension() != this->fs_alpha.dimension() || d_beta.dimension() != this->fs_beta.dimension()) {
            throw std::invalid_argument("USQOneElectronOperator::calculateFockianMatrix(OneRDM<double>, TwoRDM<double>, OneRDM<double>, TwoRDM<double>): The 2-RDM is not compatible with the one-electron operator.");
        }


        // A KISS implementation of the calculation of the Fockian matrix
        std::array<SquareMatrix<Scalar>, Components> Fs_alpha;  // alpha Fock matrices (hence the 's')
        std::array<SquareMatrix<Scalar>, Components> Fs_beta; // beta Fock matrices (hence the 's')
        for (size_t i = 0; i < Components; i++) {

            const auto& f_i_alpha = this->fs_alpha.parameters(i);  // the matrix representation of the parameters of the i-th alpha component
            const auto& f_i_beta = this->fs_beta.parameters(i); // the matrix representation of the parameters of the i-th beta component

            // Calculate the Fockian matrix for every component and add it to the array
            SquareMatrix<Scalar> F_i_alpha = SquareMatrix<Scalar>::Zero(this->fs_alpha.dimension(), this->fs_alpha.dimension());  // the alpha Fockian matrix of the i-th component
            SquareMatrix<Scalar> F_i_beta = SquareMatrix<Scalar>::Zero(this->fs_beta.dimension(), this->fs_beta.dimension());  // the beta Fockian matrix of the i-th component

            for (size_t p = 0; p < this->fs_alpha.dimension(); p++) {
                for (size_t q = 0; q < this->fs_alpha.dimension(); q++) {
                    for (size_t r = 0; r < this->fs_alpha.dimension(); r++) {
                        F_i_alpha(p,q) += f_i_alpha(q,r) * (D_alpha(p,r) + D_alpha(r,p));
                        F_i_beta(p,q) += f_i_beta(q,r) * (D_beta(p,r) + D_beta(r,p));
                    }
                }
            }  // F_i_alpha and F_i_beta elements loop. Dimensions are the same so written in 1 loop.
            Fs_alpha[i] = 0.5 * F_i_alpha;
            Fs_beta[i] = 0.5 * F_i_beta;
        }

        return Fs_alpha, Fs_beta;
    }


    /**
     *  @param D      the 1-DM (or the response 1-DM for made-variational wave function models)
     *  @param d      the 2-DM (or the response 2-DM for made-variational wave function models)
     *
     *  @return the (generalized) super-Fockian matrix
     */
    std::array<SquareRankFourTensor<Scalar>, Components> calculateSuperFockianMatrix(const OneRDM<double>& D_alpha, const TwoRDM<double>& d_alpha, const OneRDM<double>& D_beta, const TwoRDM<double>& d_beta) const {

        // Check if dimensions are compatible
        if (D_alpha.dimension() != this->fs_alpha.dimension() || D_beta.dimension() != this->fs_beta.dimension()) {
            throw std::invalid_argument("USQOneElectronOperator::calculateFockianMatrix(OneRDM<double>, TwoRDM<double>, OneRDM<double>, TwoRDM<double>): The 1-RDM is not compatible with the one-electron operator.");
        }

        if (d_alpha.dimension() != this->fs_alpha.dimension() || d_beta.dimension() != this->fs_beta.dimension()) {
            throw std::invalid_argument("USQOneElectronOperator::calculateFockianMatrix(OneRDM<double>, TwoRDM<double>, OneRDM<double>, TwoRDM<double>): The 2-RDM is not compatible with the one-electron operator.");
        }

        // A KISS implementation of the calculation of the super-Fockian matrix
        std::array<SquareRankFourTensor<Scalar>, Components> Gs_alpha;  // multiple alpha Gs, hence the 's'
        std::array<SquareRankFourTensor<Scalar>, Components> Gs_beta; // multiple beta Gs, hence the 's'
        const auto Fs_alpha, Fs_beta = this->calculateFockianMatrix(D_alpha, d_alpha, D_beta, d_beta);  // the Fockian matrices are necessary in the calculation

        for (size_t i = 0; i < Components; i++) {
            
            const auto& f_i_alpha = this->fs_alpha.parameters(i);  // the matrix representation of the parameters of the i-th alpha component
            const auto& f_i_beta = this->fs_beta.parameters(i);  // the matrix representation of the parameters of the i-th beta component

            const auto& F_i_alpha = Fs_alpha[i];  // the alpha Fockian matrix of the i-th component
            const auto& F_i_beta = Fs_beta[i];  // the beta Fockian matrix of the i-th component

            // Calculate the super-Fockian matrix for every component and add it to the array
            SquareRankFourTensor<Scalar> G_i_alpha (this->fs_alpha.dimension());
            SquareRankFourTensor<Scalar> G_i_beta (this->fs_beta.dimension());
            G_i_alpha.setZero();
            G_i_beta.setZero();
            for (size_t p = 0; p < this->fs_alpha.dimension(); p++) {
                for (size_t q = 0; q < this->fs_alpha.dimension(); q++) {
                    for (size_t r = 0; r < this->fs_alpha.dimension(); r++) {
                        for (size_t s = 0; s < this->fs_alpha.imension(); s++) {

                            if (q == r) {
                                G_i_alpha(p,q,r,s) += 2 * F_i_alpha(p,s);
                                G_i_beta(p,q,r,s) += 2 * F_i_beta(p,s);
                            }

                            G_i_alpha(p,q,r,s) -= f_i_alpha(s,p) * (D_alpha(r,q) + D_alpha(q,r));
                            G_i_beta(p,q,r,s) -= f_i_beta(s,p) * (D_beta(r,q) + D_alpha(q,r));
                        }
                    }
                }
            }  // G_i elements loop, dimensions for alpha and beta are the same and can thus be written in 1 loop.
            Gs_alpha[i] = 0.5 * G_i_alpha;
            Gs_beta[i] = 0.5 * G_i_beta;
        }

        return Gs_alpha, Gs_beta;
    }


    /**
     *  @return the dimension of the matrix representation of the parameters, i.e. the number of orbitals/sites
     */
    size_t dimension() const {
        return this->fs_alpha[0].dimension(), this->fs_beta[0].dimension();  // all the dimensions are the same, this is checked in the constructor
    }


    /**
     *  @param a        the vector
     * 
     *  @return the dot product of this second-quantized one-electron operator with the given vector
     */
    USQOneElectronOperator<Scalar, Components, Scalar, Components> dot(const Vector<Scalar, Components>& a) const {

        const auto dim = this->fs_alpha.dimension();
        const auto dim = this->fs_beta.dimension();
        USQOneElectronOperator<Scalar, Components, Scalar, Components> result_alpha, result_beta {dim};

        // Calculate the inner product
        for (size_t i = 0; i < Components; i++) {
            result_alpha += a(i) * this->fs_alpha.operator[](i);
            result_beta += a(i) * this-> fs_beta.operator[](i);           
        }

        return result_alpha, result_beta;
    }
}

} // namespace GQCP