// This file is part of GQCG-GQCP.
//
// Copyright (C) 2017-2020  the GQCG developers
//
// GQCG-GQCP is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GQCG-GQCP is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-GQCP.  If not, see <http://www.gnu.org/licenses/>.

#pragma once


#include "Mathematical/Functions/VectorSpaceArithmetic.hpp"
#include "Operator/SecondQuantized/SQOperatorStorageBase.hpp"

#include <algorithm>


namespace GQCP {


/*
 *  MARK: Implementing SQOperatorStorage
 */

/**
 *  A type that takes care of storing the matrix elements of one- and two-electron operators and their components, which adds functionality for vector space arithmetic out of the box.
 * 
 *  @tparam _MatrixRepresentation       The type used to represent the set of parameters/matrix elements/integrals for one component of a second-quantized operator.
 *  @tparam _Vectorizer                 The type of the vectorizer that relates a one-dimensional storage of matrix representations to the tensor structure of the second-quantized operators. This allows for a distinction between scalar operators (such as the kinetic energy or Coulomb operator), vector operators (such as the spin operator) and matrix/tensor operators (such as quadrupole and multipole operators).
 *  @tparam _DerivedOperator            The type of the operator that derives from this class, enabling CRTP and compile-time polymorphism.
 */
template <typename _MatrixRepresentation, typename _Vectorizer, typename _DerivedOperator>
class SQOperatorStorage:
    public SQOperatorStorageBase<_MatrixRepresentation, _Vectorizer, _DerivedOperator>,
    public VectorSpaceArithmetic<typename OperatorTraits<_DerivedOperator>::DerivedOperator, typename OperatorTraits<_DerivedOperator>::Scalar> {
public:
    // The type used to represent the set of parameters/matrix elements/integrals for one component of a second-quantized operator.
    using MatrixRepresentation = _MatrixRepresentation;

    // The type of the operator that derives from this class, enabling CRTP and compile-time polymorphism.
    using DerivedOperator = _DerivedOperator;

    // The scalar type used for a single parameter/matrix element/integral: real or complex.
    using Scalar = typename OperatorTraits<DerivedOperator>::Scalar;

    // The type of the vectorizer that relates a one-dimensional storage of matrix representations to the tensor structure of the second-quantized operators. This allows for a distinction between scalar operators (such as the kinetic energy or Coulomb operator), vector operators (such as the spin operator) and matrix/tensor operators (such as quadrupole and multipole operators).
    using Vectorizer = _Vectorizer;

    // The type of the final derived operator, enabling CRTP and compile-time polymorphism. VectorSpaceArithmetic (and other functionality) should be implemented on the **final** deriving class, not on intermediate classes. In the current design of the SQOperator classes, there is just one intermediate class, so finding out the final derived class is easy.
    using FinalOperator = typename OperatorTraits<DerivedOperator>::DerivedOperator;

    // The type that corresponds to the scalar version of the final second-quantized operator operator type.
    using ScalarFinalOperator = typename OperatorTraits<FinalOperator>::ScalarOperator;


public:
    /*
     *  MARK: Constructors
     */

    // Inherit `SQOperatorStorageBase`'s constructors.
    using SQOperatorStorageBase<MatrixRepresentation, Vectorizer, DerivedOperator>::SQOperatorStorageBase;


    /*
     *  MARK: Conforming to `VectorSpaceArithmetic`
     */

    /**
     *  Addition-assignment.
     */
    FinalOperator& operator+=(const FinalOperator& rhs) override {

        // Use the STL to implement element-wise addition.
        std::transform(this->array.elements().begin(), this->array.elements().end(),
                       rhs.array.elements().begin(), this->array.elements().begin(),
                       [](const MatrixRepresentation& M_lhs, const MatrixRepresentation& M_rhs) {
                           return M_lhs.Eigen() + M_rhs.Eigen();
                       });

        return static_cast<FinalOperator&>(*this);
    }


    /**
     *  Scalar multiplication-assignment.
     */
    FinalOperator& operator*=(const Scalar& a) override {

        // Use the STL to implement element-wise scalar multiplication.
        std::transform(this->array.elements().begin(), this->array.elements().end(),
                       this->array.elements().begin(), [a](const MatrixRepresentation& M) { return M * a; });

        return static_cast<FinalOperator&>(*this);
    }


    /*
     *  MARK: Other arithmetic
     */

    /**
     *  Calculate the scalar one-electron operator that is the result of the vector dot product between this operator and the given vector. Note that this method is only available for vector-like one-electron operators.
     * 
     *  @param v        The vector with which this operator should be dot-multiplied.
     * 
     *  @return The vector dot product between this operator and the given vector.
     */
    template <typename Z = Vectorizer>
    ScalarFinalOperator dot(const VectorX<Scalar>& v) const {

        if (this->array.vectorizer().numberOfElements() != v.size()) {
            throw std::invalid_argument("SimpleSQOneElectronOperator(const StorageArray<Scalar, Vectorizer>&): The dimension of the given vector is incompatible with this vector operator.)");
        }

        const auto dimension = this->numberOfOrbitals();
        ScalarFinalOperator result = ScalarFinalOperator::Zero(dimension);  // Initializes a scalar-like operator with parameters that are zero.

        // Calculate the dot/inner product of two vectors.
        for (size_t i = 0; i < this->numberOfComponents(); i++) {
            result += v(i) * (*this)(i);
        }

        return result;
    }
};


}  // namespace GQCP
