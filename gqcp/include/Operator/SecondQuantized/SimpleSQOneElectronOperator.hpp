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


#include "Basis/Transformations/BasisTransformable.hpp"
#include "Basis/Transformations/JacobiRotatable.hpp"
#include "DensityMatrix/OneDM.hpp"
#include "Mathematical/Functions/VectorSpaceArithmetic.hpp"
#include "Mathematical/Representation/SquareMatrix.hpp"
#include "Mathematical/Representation/StorageArray.hpp"
#include "Utilities/CRTP.hpp"
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


/**
 *  A second-quantized one-electron operator whose parameters are described by a single matrix.
 * 
 *  This class is used as a base class for `RSQOneElectronOperator` and `GSQOneElectronOperator`, since they both admit parameter representations using a single matrix, as opposed to `USQOneElectronOperator`, which uses separate alpha- and beta- matrices. The word 'simple' is used here as an antonym for 'compound'.
 * 
 *  @tparam _Scalar                 The scalar type used for a single parameter: real or complex.
 *  @tparam _Vectorizer             The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of one-electron operators. This allows for a distinction between scalar operators (such as the kinetic energy operator), vector operators (such as the spin operator) and matrix/tensor operators (such as quadrupole and multipole operators).
 *  @tparam _DerivedOperator        The type of the operator that derives from this class, enabling CRTP and compile-time polymorphism.
 */
template <typename _Scalar, typename _Vectorizer, typename _DerivedOperator>
class SimpleSQOneElectronOperator:
    public VectorSpaceArithmetic<_DerivedOperator, _Scalar>,
    public BasisTransformable<_DerivedOperator, typename OperatorTraits<_DerivedOperator>::TM>,
    public JacobiRotatable<_DerivedOperator> {
public:
    // The scalar type used for a single parameter: real or complex.
    using Scalar = _Scalar;

    // The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of one-electron operators. This allows for a distinction between scalar operators (such as the kinetic energy operator), vector operators (such as the spin operator) and matrix/tensor operators (such as quadrupole and multipole operators).
    using Vectorizer = _Vectorizer;

    // The type of the operator that derives from this class, enabling CRTP and compile-time polymorphism.
    using DerivedOperator = _DerivedOperator;

    // The matrix representation of the parameters of (one of the components of) the one-electron operator.
    using MatrixRepresentation = SquareMatrix<Scalar>;

    // The type of 'this'.
    using Self = SimpleSQOneElectronOperator<Scalar, Vectorizer, DerivedOperator>;

    // The type that corresponds to the scalar version of the derived one-electron operator type.
    using ScalarDerivedOperator = typename OperatorTraits<DerivedOperator>::ScalarOperator;

    // The type of transformation matrix that is naturally associated to the derived one-electron operator.
    using TM = typename OperatorTraits<DerivedOperator>::TM;

    // The type of density matrix that is naturally associated to the derived one-electron operator.
    using Derived1DM = typename OperatorTraits<DerivedOperator>::OneDM;


private:
    // The array that takes care of the storage of the operator's matrix elements.
    StorageArray<MatrixRepresentation, Vectorizer> array;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  Construct a one-electron operator from a(n) (storage) array.
     * 
     *  @param array                A storage array that contains the matrix representations of all the components of this operator.
     */
    SimpleSQOneElectronOperator(const StorageArray<MatrixRepresentation, Vectorizer>& array) :
        array {array} {

        const auto first_dimension = array.elements()[0].dimension();

        for (const auto& parameters : array.elements()) {
            if (parameters.dimension() != first_dimension) {
                throw std::invalid_argument("SimpleSQOneElectronOperator(const StorageArray<MatrixRepresentation, Vectorizer>& array): The dimensions of the matrix representations must be equal.");
            }
        }
    }


    /**
     *  Construct a one-electron operator from one matrix representation.
     * 
     *  @param parameters           The matrix representation of the one-electron integrals/parameters.
     */
    template <typename Z = Vectorizer>
    SimpleSQOneElectronOperator(const MatrixRepresentation& parameters, typename std::enable_if<std::is_same<Z, ScalarVectorizer>::value>::type* = 0) :
        SimpleSQOneElectronOperator(StorageArray<MatrixRepresentation, ScalarVectorizer>({parameters}, ScalarVectorizer())) {}


    /**
     *  Construct a one-electron operator from a set of three matrix representations, for each of the operator's components.
     * 
     *  @param parameters           A set of three matrix representations of the one-electron integrals/parameters, one for each of the operator's components.
     */
    template <typename Z = Vectorizer>
    SimpleSQOneElectronOperator(const std::vector<MatrixRepresentation>& parameters, typename std::enable_if<std::is_same<Z, VectorVectorizer>::value>::type* = 0) :
        SimpleSQOneElectronOperator(StorageArray<MatrixRepresentation, VectorVectorizer>(parameters, VectorVectorizer({3}))) {}


    /**
     *  Construct a one-electron operator with parameters that are zero.
     * 
     *  @param dim          The dimension of the matrix representation of the parameters, i.e. the number of orbitals/sites.
     */
    SimpleSQOneElectronOperator(const size_t dim, const Vectorizer& vectorizer) :
        SimpleSQOneElectronOperator(StorageArray<MatrixRepresentation, Vectorizer> {MatrixRepresentation::Zero(dim), vectorizer}) {}


    /**
     *  Construct a one-electron operator with parameters that are zero.
     * 
     *  @param dim          The dimension of the matrix representation of the parameters, i.e. the number of orbitals/sites.
     */
    template <typename Z = Vectorizer>
    SimpleSQOneElectronOperator(const size_t dim, typename std::enable_if<std::is_same<Z, ScalarVectorizer>::value>::type* = 0) :
        SimpleSQOneElectronOperator(dim, ScalarVectorizer()) {}


    /**
     *  Construct a one-electron operator with parameters that are zero.
     * 
     *  @param dim          The dimension of the matrix representation of the parameters, i.e. the number of orbitals/sites.
     */
    template <typename Z = Vectorizer>
    SimpleSQOneElectronOperator(const size_t dim, typename std::enable_if<std::is_same<Z, VectorVectorizer>::value>::type* = 0) :
        SimpleSQOneElectronOperator(dim, VectorVectorizer({3})) {}


    /**
     *  The default constructor.
     */
    SimpleSQOneElectronOperator() :
        SimpleSQOneElectronOperator(0, Vectorizer()) {}


    /*
     *  MARK: Parameter access
     */

    /**
     *  @return A vector of read-only matrix representions of the parameters/integrals of this operator.
     */
    const std::vector<MatrixRepresentation>& allParameters() const { return this->array.elements(); }

    /**
     *  @return A vector of writable matrix representions of the parameters/integrals of this operator.
     */
    std::vector<MatrixRepresentation>& allParameters() { return this->array.elements(); }

    /**
     *  @param indices      A set of coordinates that accesses this one-electron operator.
     * 
     *  @return A read-only matrix representation of the parameters/integrals of one of the tensor components of this operator.
     */
    template <typename... Indices>
    const MatrixRepresentation& parameters(const Indices&... indices) const { return this->array(indices...); }


    /**
     *  @param indices      A set of coordinates that accesses this one-electron operator.
     * 
     *  @return A writable matrix representation of the parameters/integrals of one of the tensor components of this operator.
     */
    template <typename... Indices>
    MatrixRepresentation& parameters(const Indices&... indices) { return this->array(indices...); }


    /**
     *  @param indices      A set of coordinates that accesses this one-electron operator.
     * 
     *  @return The component of this one-electron operator that corresponds to the given coordinate indices.
     */
    template <typename... Indices>
    ScalarDerivedOperator operator()(const Indices&... indices) const {

        // Access the underlying array's storage, and wrap the result in a scalar form of the derived second-quantized operator. The correct scalar derived operator can be instantiated since the scalar derived operator type has a special constructor that expects just one matrix.
        return ScalarDerivedOperator {this->array(indices...)};
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
     *  MARK: Calculations
     */

    /**
     *  Calculate the expectation value of this one-electron operator.
     * 
     *  @param D                The 1-DM that represents the wave function.
     *
     *  @return The expectation value of all components of the one-electron operator.
     */
    StorageArray<Scalar, Vectorizer> calculateExpectationValue(const Derived1DM& D) const {

        if (this->numberOfOrbitals() != D.numberOfOrbitals()) {
            throw std::invalid_argument("SimpleSQOneElectronOperator::calculateExpectationValue(const OneDM<Scalar>&): The given 1-DM is not compatible with the one-electron operator.");
        }

        // Calculate the expectation value for every component of the operator.
        const auto& parameters = this->allParameters();
        std::vector<Scalar> expectation_values(this->numberOfComponents());  // zero-initialize the vector with a number of elements
        for (size_t i = 0; i < this->numberOfComponents(); i++) {
            expectation_values[i] = (parameters[i] * D).trace();
        }

        return StorageArray<Scalar, Vectorizer> {expectation_values, this->array.vectorizer()};
    }


    /*
     *  MARK: Conforming to VectorSpaceArithmetic
     */

    /**
     *  Addition-assignment.
     */
    DerivedOperator& operator+=(const DerivedOperator& rhs) override {

        // Use the STL to implement element-wise addition.
        std::transform(this->array.elements().begin(), this->array.elements().end(),
                       rhs.array.elements().begin(), this->array.elements().begin(),
                       std::plus<MatrixRepresentation>());

        return static_cast<DerivedOperator&>(*this);
    }


    /**
     *  Scalar multiplication-assignment.
     */
    DerivedOperator& operator*=(const Scalar& a) override {

        // Use the STL to implement element-wise scalar multiplication.
        std::transform(this->array.elements().begin(), this->array.elements().end(),
                       this->array.elements().begin(), [a](const MatrixRepresentation& M) { return M * a; });

        return static_cast<DerivedOperator&>(*this);
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
    ScalarDerivedOperator dot(const VectorX<Scalar>& v) const {

        if (this->array.vectorizer().numberOfElements() != v.size()) {
            throw std::invalid_argument("SimpleSQOneElectronOperator(const StorageArray<Scalar, Vectorizer>&): The dimension of the given vector is incompatible with this vector operator.)");
        }

        const auto dimension = this->numberOfOrbitals();
        ScalarDerivedOperator result {dimension};  // initializes a scalar one-electron operator with parameters that are zero

        // Calculate the dot/inner product of two vectors.
        for (size_t i = 0; i < this->numberOfComponents(); i++) {
            result += v(i) * (*this)(i);
        }

        return result;
    }


    /*
     *  MARK: Conforming to BasisTransformable
     */

    /**
     *  Apply the basis transformation and return the resulting one-electron integrals.
     * 
     *  @param transformation_matrix        The type that encapsulates the basis transformation coefficients.
     * 
     *  @return The basis-transformed one-electron integrals.
     */
    DerivedOperator transformed(const TM& transformation_matrix) const override {

        // Calculate the basis transformation for every component of the operator.
        const auto& parameters = this->allParameters();
        auto result = this->allParameters();

        for (size_t i = 0; i < this->numberOfComponents(); i++) {
            result[i] = transformation_matrix.adjoint() * (parameters[i]) * transformation_matrix;
        }

        return DerivedOperator {StorageArray<MatrixRepresentation, Vectorizer>(result, this->array.vectorizer())};
    }


    // Allow the `rotate` method from `BasisTransformable`, since there's also a `rotate` from `JacobiRotatable`.
    using BasisTransformable<DerivedOperator, TM>::rotate;

    // Allow the `rotated` method from `BasisTransformable`, since there's also a `rotated` from `JacobiRotatable`.
    using BasisTransformable<DerivedOperator, TM>::rotated;


    /*
     *  MARK: Conforming to JacobiRotatable
     */

    /**
     *  Apply the Jacobi rotation and return the result.
     * 
     *  @param jacobi_parameters        The Jacobi rotation parameters.
     * 
     *  @return The jacobi-transformed object.
     */
    DerivedOperator rotated(const JacobiRotationParameters& jacobi_parameters) const override {

        // Use Eigen's Jacobi module to apply the Jacobi rotations directly (cfr. T.adjoint() * M * T).
        const auto p = jacobi_parameters.p();
        const auto q = jacobi_parameters.q();
        const auto jacobi_rotation = jacobi_parameters.Eigen();

        // Calculate the basis transformation for every component of the operator.
        auto result = this->allParameters();
        for (size_t i = 0; i < this->numberOfComponents(); i++) {
            result[i].applyOnTheLeft(p, q, jacobi_rotation.adjoint());
            result[i].applyOnTheRight(p, q, jacobi_rotation);
        }

        return DerivedOperator {StorageArray<MatrixRepresentation, Vectorizer>(result, this->array.vectorizer())};
    }

    // Allow the `rotate` method from `JacobiRotatable`, since there's also a `rotate` from `BasisTransformable`.
    using JacobiRotatable<DerivedOperator>::rotate;
};


}  // namespace GQCP
