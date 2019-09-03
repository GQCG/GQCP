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


#include "Basis/TransformationMatrix.hpp"
#include "Mathematical/Representation/ChemicalMatrix.hpp"
#include "Mathematical/ScalarFunction.hpp"
#include "OrbitalOptimization/JacobiRotationParameters.hpp"
#include "RDM/OneRDM.hpp"
#include "RDM/TwoRDM.hpp"
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
class SQOneElectronOperator {
public:

    using Scalar = _Scalar;
    static constexpr auto Components = _Components;


private:
    std::array<ChemicalMatrix<Scalar>, Components> fs;  // all the matrix representations (hence the s) of the parameters (integrals) of the different components of this second-quantized operator


public:

    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param fs           all the matrix representations (hence the s) of the parameters (integrals) of the different components of this second-quantized operator
     */
    SQOneElectronOperator(const std::array<ChemicalMatrix<Scalar>, Components>& fs) : 
        fs (fs)
    {
        // Check if the given matrix representations have the same dimensions
        const auto dimension_of_first = this->fs[0].dimension();
        for (size_t i = 1; i < Components; i++) {

            const auto dimension_of_ith = this->fs[i].dimension();
            if (dimension_of_first != dimension_of_ith) {
                throw std::invalid_argument("SQOneElectronOperator(const std::array<ChemicalMatrix<Scalar>, Components>&): The given matrix representations do not have the same dimensions.");
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
            this->fs[i] = ChemicalMatrix<Scalar>::Zero(dim, dim);
        }
    }


    /**
     *  Default constructor: construct a one-electron operator with parameters that are zero
     */
    SQOneElectronOperator() :
        SQOneElectronOperator(0)  // dimensions of the representations are zero
    {}


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return the dimension of the matrix representation of the parameters, i.e. the number of orbitals/sites
     */
    size_t dimension() const {
        return this->fs[0].dimension();  // all the dimensions are the same, this is checked in the constructor
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
    const std::array<ChemicalMatrix<Scalar>, Components>& allParameters() const {
        return this->fs;
    }


    /**
     *  @return writable matrix representations of all the parameters (integrals) of the different components of this second-quantized operator
     */
    std::array<ChemicalMatrix<Scalar>, Components>& allParameters() {
        return this->fs;
    }


    /**
     *  @param i            the index of the component
     * 
     *  @return a read-only the matrix representation of the parameters (integrals) of one of the the different components of this second-quantized operator
     */
    const ChemicalMatrix<Scalar>& parameters(const size_t i = 0) const {
        return this->fs[i];
    }


    /**
     *  @param i            the index of the component
     * 
     *  @return a writable matrix representation of the parameters (integrals) of one of the the different components of this second-quantized operator
     */
    ChemicalMatrix<Scalar>& parameters(const size_t i = 0) {
        return this->fs[i];
    }


    /**
     *  In-place transform the operator to another basis
     * 
     *  @param T                            the transformation matrix
     */
    void transform(const TransformationMatrix<Scalar>& T) {

        // Transform the matrix representations of the components
        for (auto& f : this->allParameters()) {
            f.basisTransformInPlace(T);
        }
    }


    /**
     *  In-place rotate the operator to another basis
     * 
     *  @param U                            the (unitary) rotation matrix
     */
    void rotate(const TransformationMatrix<Scalar>& U) {

        // Transform the matrix representations of the components
        for (auto& f : this->allParameters()) {
            f.basisRotateInPlace(U);
        }
    }


    /**
     *  In-place rotate the operator using a unitary Jacobi rotation matrix constructed from the Jacobi rotation parameters
     * 
     *  @param jacobi_rotation_parameters       the Jacobi rotation parameters (p, q, angle) that are used to specify a Jacobi rotation: we use the (cos, sin, -sin, cos) definition for the Jacobi rotation matrix. See transform() for how the transformation matrix between the two bases should be represented
     */
    void rotate(const JacobiRotationParameters& jacobi_rotation_parameters) {

        // Transform the matrix representations of the components
        for (auto& f : this->allParameters()) {
            f.basisRotateInPlace(jacobi_rotation_parameters);
        }
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
        std::array<ChemicalMatrix<typename Z::Valued>, Components> F_evaluated;  // components are not initialized here

        // Evaluate all components at the given x
        for (size_t i = 0; i < Components; i++) {
            F_evaluated[i] = ChemicalMatrix<typename Z::Valued>::Zero(this->dimension(), this->dimension());  // initialize to zero

            for (size_t m = 0; m < this->dimension(); m++) {
                for (size_t n = 0; n < this->dimension(); n++) {
                    const auto F_i_mn = this->parameters(i)(m,n);  // (m,n)-th element of the i-th component
                    F_evaluated[i](m,n) = F_i_mn.operator()(x);  // evaluate the ScalarFunction
                }
            }
        }

        return SQOneElectronOperator<typename Z::Valued, Components>(F_evaluated);
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
            throw std::invalid_argument("SQOneElectronOperator::calculateFockianMatrix(OneRDM<double>, TwoRDM<double>): The 1-RDM is not compatible with the one-electron operator.");
        }

        if (d.dimension() != this->dimension()) {
            throw std::invalid_argument("SQOneElectronOperator::calculateFockianMatrix(OneRDM<double>, TwoRDM<double>): The 2-RDM is not compatible with the one-electron operator.");
        }


        // A KISS implementation of the calculation of the Fockian matrix
        std::array<SquareMatrix<Scalar>, Components> Fs;  // Fock matrices (hence the 's')
        for (size_t i = 0; i < Components; i++) {

            const auto& f_i = this->parameters(i);  // the matrix representation of the parameters of the i-th component

            // Calculate the Fockian matrix for every component and add it to the array
            SquareMatrix<Scalar> F_i = SquareMatrix<Scalar>::Zero(this->dimension(), this->dimension());  // the Fockian matrix of the i-th component
            for (size_t p = 0; p < this->dimension(); p++) {
                for (size_t q = 0; q < this->dimension(); q++) {

                    for (size_t r = 0; r < this->dimension(); r++) {
                        F_i(p,q) += f_i(q,r) * (D(p,r) + D(r,p));
                    }

                }
            }  // F_i elements loop
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
            
            const auto& f_i = this->parameters(i);  // the matrix representation of the parameters of the i-th component
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

                            G_i(p,q,r,s) -= f_i(s,p) * (D(r,q) + D(q,r));

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

    auto fs = op.allParameters();
    for (auto& f : fs) {
        f *= scalar;
    }

    return SQOneElectronOperator<ResultScalar, Components>(fs);
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
SQOneElectronOperator<Scalar, Components> operator-(const SQOneElectronOperator<Scalar, Components>& op) {

    return (-1.0) * op;  // negation is scalar multiplication with (-1.0)
}


/**
 *  Subtract two one-electron operators by adding their parameters
 * 
 *  @tparam LHSScalar           the scalar type of the left-hand side
 *  @tparam RHSScalar           the scalar type of the right-hand side
 *  @tparam Components          the number of components of the one-electron operators
 * 
 *  @param lhs                  the left-hand side
 *  @param rhs                  the right-hand side
 */
template <typename LHSScalar, typename RHSScalar, size_t Components>
auto operator-(const SQOneElectronOperator<LHSScalar, Components>& lhs, const SQOneElectronOperator<RHSScalar, Components>& rhs) -> SQOneElectronOperator<sum_t<LHSScalar, RHSScalar>, Components> {

    return lhs + (-rhs);
}


}  // namespace GQCP
