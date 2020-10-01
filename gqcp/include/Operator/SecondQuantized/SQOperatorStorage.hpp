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
#include "Mathematical/Representation/StorageArray.hpp"
#include "Utilities/type_traits.hpp"

#include <algorithm>


namespace GQCP {


/*
 *  MARK: Operator traits
 */

/**
 *  A type that provides compile-time information on operators that is otherwise not accessible through a public class alias.
 */
template <typename Operator>
class OperatorTraits {};


/*
 *  MARK: Implementing SQOperatorStorage
 */

/**
 *  A type that takes care of storing the matrix elements of one- and two-electron operators and their components.
 * 
 *  @tparam _MatrixRepresentation       The type used to represent the set of parameters/matrix elements/integrals for one component of a second-quantized operator.
 *  @tparam _Vectorizer                 The type of the vectorizer that relates a one-dimensional storage of matrix representations to the tensor structure of the second-quantized operators. This allows for a distinction between scalar operators (such as the kinetic energy or Coulomb operator), vector operators (such as the spin operator) and matrix/tensor operators (such as quadrupole and multipole operators).
 *  @tparam _DerivedOperator            The type of the operator that derives from this class, enabling CRTP and compile-time polymorphism.
 */
template <typename _MatrixRepresentation, typename _Vectorizer, typename _DerivedOperator>
class SQOperatorStorage:
    public VectorSpaceArithmetic<typename OperatorTraits<_DerivedOperator>::DerivedOperator, typename OperatorTraits<_DerivedOperator>::Scalar> {
public:
    // The type used to represent the set of parameters/matrix elements/integrals for one component of a second-quantized operator.
    using MatrixRepresentation = _MatrixRepresentation;

    // The scalar type used for a single parameter/matrix element/integral: real or complex.
    using Scalar = typename MatrixRepresentation::Scalar;

    // The type of the vectorizer that relates a one-dimensional storage of matrix representations to the tensor structure of the second-quantized operators. This allows for a distinction between scalar operators (such as the kinetic energy or Coulomb operator), vector operators (such as the spin operator) and matrix/tensor operators (such as quadrupole and multipole operators).
    using Vectorizer = _Vectorizer;

    // The type of the operator that derives from this class, enabling CRTP and compile-time polymorphism.
    using DerivedOperator = _DerivedOperator;

    // The type of the final derived operator, enabling CRTP and compile-time polymorphism. VectorSpaceArithmetic (and other functionality) should be implemented on the **final** deriving class, not on intermediate classes. In the current design of the SQOperator classes, there is just one intermediate class, so finding out the final derived class is easy.
    using FinalOperator = typename OperatorTraits<_DerivedOperator>::DerivedOperator;

    // The type that corresponds to the scalar version of the final second-quantized operator operator type.
    using ScalarFinalOperator = typename OperatorTraits<FinalOperator>::ScalarOperator;


protected:
    // The array that takes care of the storage of the second quantized operator's matrix elements.
    StorageArray<MatrixRepresentation, Vectorizer> array;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  Construct a second-quantized operator storage from a storage array.
     * 
     *  @param array                A storage array that contains the matrix representations of all the components of the operator.
     */
    SQOperatorStorage(const StorageArray<MatrixRepresentation, Vectorizer>& array) :
        array {array} {

        const auto first_dimension = array.elements()[0].dimension();

        for (const auto& parameters : array.elements()) {
            if (parameters.dimension() != first_dimension) {
                throw std::invalid_argument("SQOperatorStorage(const StorageArray<MatrixRepresentation, Vectorizer>& array): The dimensions of the matrix representations must be equal.");
            }
        }
    }


    /**
     *  Construct a second-quantized operator storage from one matrix representation.
     * 
     *  @param parameters           The matrix representation of operator's parameters/matrix elements/integrals.
     */
    template <typename Z = Vectorizer>
    SQOperatorStorage(const MatrixRepresentation& parameters, typename std::enable_if<std::is_same<Z, ScalarVectorizer>::value>::type* = 0) :
        SQOperatorStorage(StorageArray<MatrixRepresentation, ScalarVectorizer>({parameters}, ScalarVectorizer())) {}


    /**
     *  Construct a second-quantized operator storage from a set of three matrix representations, for each of the operator's components.
     * 
     *  @param parameters           A set of three matrix representations of the one-electron parameters/matrix elements/integrals, one for each of the operator's components.
     */
    template <typename Z = Vectorizer>
    SQOperatorStorage(const std::vector<MatrixRepresentation>& parameters, typename std::enable_if<std::is_same<Z, VectorVectorizer>::value>::type* = 0) :
        SQOperatorStorage(StorageArray<MatrixRepresentation, VectorVectorizer>(parameters, VectorVectorizer({3}))) {}


    /**
     *  Construct a second-quantized operator storage with parameters/matrix elements/integrals that are zero.
     * 
     *  @param dim          The dimension of the matrix representation of the parameters/matrix elements/integrals, i.e. the number of orbitals/sites.
     */
    SQOperatorStorage(const size_t dim, const Vectorizer& vectorizer) :
        SQOperatorStorage(StorageArray<MatrixRepresentation, Vectorizer> {MatrixRepresentation::Zero(dim), vectorizer}) {}


    /**
     *  Construct a scalar second-quantized operator storage with parameters/matrix elements/integrals that are zero.
     * 
     *  @param dim          The dimension of the matrix representation of the parameters/matrix elements/integrals, i.e. the number of orbitals/sites.
     */
    template <typename Z = Vectorizer>
    SQOperatorStorage(const size_t dim, typename std::enable_if<std::is_same<Z, ScalarVectorizer>::value>::type* = 0) :
        SQOperatorStorage(dim, ScalarVectorizer()) {}


    /**
     *  Construct a vector second-quantized operator storage with parameters/matrix elements/integrals that are zero.
     * 
     *  @param dim          The dimension of the matrix representation of the parameters/matrix elements/integrals, i.e. the number of orbitals/sites.
     */
    template <typename Z = Vectorizer>
    SQOperatorStorage(const size_t dim, typename std::enable_if<std::is_same<Z, VectorVectorizer>::value>::type* = 0) :
        SQOperatorStorage(dim, VectorVectorizer({3})) {}


    /**
     *  The default constructor.
     */
    SQOperatorStorage() :
        SQOperatorStorage(0, Vectorizer()) {}


    /*
     *  MARK: Parameter access
     */

    /**
     *  @return A vector of read-only matrix representions of the parameters/matrix elements/integrals of this operator.
     */
    const std::vector<MatrixRepresentation>& allParameters() const { return this->array.elements(); }

    /**
     *  @return A vector of writable matrix representions of the parameters/matrix elements/integrals of this operator.
     */
    std::vector<MatrixRepresentation>& allParameters() { return this->array.elements(); }

    /**
     *  @param indices      A set of coordinates that accesses this one-electron operator.
     * 
     *  @return A read-only matrix representation of the parameters/matrix elements/integrals of one of the tensor components of this operator.
     */
    template <typename... Indices>
    const MatrixRepresentation& parameters(const Indices&... indices) const { return this->array(indices...); }


    /**
     *  @param indices      A set of coordinates that accesses this one-electron operator.
     * 
     *  @return A writable matrix representation of the parameters/matrix elements/integrals of one of the tensor components of this operator.
     */
    template <typename... Indices>
    MatrixRepresentation& parameters(const Indices&... indices) { return this->array(indices...); }


    /**
     *  @param indices      A set of coordinates that accesses this one-electron operator.
     * 
     *  @return The component of this one-electron operator that corresponds to the given coordinate indices.
     */
    template <typename... Indices>
    ScalarFinalOperator operator()(const Indices&... indices) const {

        // Access the underlying array's storage, and wrap the result in a scalar form of the derived second-quantized operator. The correct scalar derived operator can be instantiated since the scalar derived operator type has a special constructor that expects just one matrix.
        return ScalarFinalOperator {this->array(indices...)};
    }


    /*
     *  MARK: General info
     */

    /**
     *  @return The number of components of this one-electron operator. For scalar operators (such as the kinetic energy), this is 1, for vector operators (like the dipole operator), this is 3, etc.
     */
    size_t numberOfComponents() const { return this->array.vectorizer().numberOfElements(); }

    /**
     *  @return The number of orbitals this one-electron operator is quantized in. For 'restricted' operators, this is the number of spatial orbitals, for 'general' operators, this is the number of spinors.
     */
    size_t numberOfOrbitals() const { return this->array.elements()[0].dimension(); /* all the dimensions are the same, this is checked in the constructor */ }


    /*
     *  MARK: Conforming to VectorSpaceArithmetic
     */

    /**
     *  Addition-assignment.
     */
    FinalOperator& operator+=(const FinalOperator& rhs) override {

        // Use the STL to implement element-wise addition.
        std::transform(this->array.elements().begin(), this->array.elements().end(),
                       rhs.array.elements().begin(), this->array.elements().begin(),
                       std::plus<MatrixRepresentation>());

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
        ScalarFinalOperator result {dimension};  // initializes a scalar one-electron operator with parameters that are zero

        // Calculate the dot/inner product of two vectors.
        for (size_t i = 0; i < this->numberOfComponents(); i++) {
            result += v(i) * (*this)(i);
        }

        return result;
    }
};


}  // namespace GQCP
