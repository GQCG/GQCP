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
template <typename _Scalar, size_t _Components>
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
    size_t get_dim_alpha() const { return this->fs_alpha[0].dimension(); }
    size_t get_dim_beta() const { return this->fs_beta[0].dimension(); }
    size_t get_K_alpha() const { return this->fs_alpha[0].dimension(); }
    size_t get_K_beta() const { return this->fs_beta[0].dimension(); }


    /*
     *  OPERATORS
     */
    /**
     *  @param i            the index
     * 
     *  @return the i-th component of this operator, for both spin components
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

        if (this->fs_alpha[0].dimension() != D_alpha.dimension() || this->fs_beta[0].dimension() != D_beta.dimension()) {
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
    USQOneElectronOperator<Scalar, Components> dot(const Vector<Scalar, Components>& a) const {

        const auto dim = this->fs_alpha[0].dimension();
        const auto dim = this->fs_beta[0].dimension();
        USQOneElectronOperator<Scalar, Components> result_alpha, result_beta {dim};

        // Calculate the inner product
        for (size_t i = 0; i < Components; i++) {
            result_alpha += a(i) * this->fs_alpha.operator[](i);
            result_beta += a(i) * this-> fs_beta.operator[](i);           
        }

        return result_alpha, result_beta;
    }


    /**
     *  @param i            the index of the component
     * 
     *  @return a read-only the matrix representation of the parameters (integrals) of one of the the different components of this second-quantized operator
     */
    const QCMatrix<Scalar>& parameters(const size_t i = 0) const {
        return this->fs_alpha[i], this->fs_beta[i];
    }


    /**
     *  @param i            the index of the component
     * 
     *  @return a writable matrix representation of the parameters (integrals) of one of the the different components of this second-quantized operator
     */
    QCMatrix<Scalar>& parameters(const size_t i = 0) {
        return this->fs_alpha[i], this->fs_beta[i];
    }


    /**
     *  In-place rotate the operator to another basis
     * 
     *  @param U                            the (unitary) rotation matrix
     */
    void rotate(const TransformationMatrix<Scalar>& U) {

        // Transform the matrix representations of the components
        for (auto& f_alpha : this->fs_alpha.allParameters()) {
            f_alpha.basisRotateInPlace(U);
        }
        for (auto& f_beta : this->fs_beta.allParameters())) {
            f_beta.basisRotateInPlace(U);
        }
    }


    /**
     *  In-place rotate the operator using a unitary Jacobi rotation matrix constructed from the Jacobi rotation parameters
     * 
     *  @param jacobi_rotation_parameters       the Jacobi rotation parameters (p, q, angle) that are used to specify a Jacobi rotation: we use the (cos, sin, -sin, cos) definition for the Jacobi rotation matrix
     */
    void rotate(const JacobiRotationParameters& jacobi_rotation_parameters) {

        // Transform the matrix representations of the components
        for (auto& f_alpha : this->fs_alpha.allParameters()) {
            f_alpha.basisRotateInPlace(jacobi_rotation_parameters);
        }
        for (auto& f_beta : this->fs_beta.allParameters()) {
            f_beta.basisRotateInPlace(jacobi_rotation_parameters);
        }
    }


    /**
     *  In-place transform the operator to another basis
     * 
     *  @param T                            the transformation matrix
     */
    void transform(const TransformationMatrix<Scalar>& T) {

        // Transform the matrix representations of the components
        for (auto& f_alpha : this->fs_alpha.allParameters()) {
            f_alpha.basisTransformInPlace(T);
        }
        for (auto& f_beta : this->fs_beta.allParameters()) {
            f_beta.basisTransformInPlace(T);
        }
    }    
};

} // namespace GQCP