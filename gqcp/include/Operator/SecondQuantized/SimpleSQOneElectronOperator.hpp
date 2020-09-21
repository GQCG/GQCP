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


#include "Mathematical/Representation/QCMatrix.hpp"
#include "Mathematical/Representation/StorageArray.hpp"
#include "Processing/DensityMatrices/OneDM.hpp"
#include "Utilities/type_traits.hpp"


namespace GQCP {


/**
 *  @brief A second-quantized one-electron operator whose parameters are described by a single matrix.
 * 
 *  This class is used as a base class for `RSQOneElectronOperator` and `GSQOneElectronOperator`, since they both admit parameter representations using a single matrix, as opposed to `USQOneElectronOperator`, which uses separate alpha- and beta- matrices.
 * 
 *  @tparam _Scalar         The scalar type used for a single parameter: real or complex.
 *  @tparam _Vectorizer     The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of one-electron operators. This allows for a distinction between scalar operators (such as the kinetic energy operator), vector operators (such as the spin operator) and matrix/tensor operators (such as quadrupole and multipole operators).
 */
template <typename _Scalar, typename _Vectorizer>
class SimpleSQOneElectronOperator {
public:
    // The scalar type used for a single parameter: real or complex.
    using Scalar = _Scalar;

    // The type of the vectorizer that relates a one-dimensional storage of matrices to the tensor structure of one-electron operators. This allows for a distinction between scalar operators (such as the kinetic energy operator), vector operators (such as the spin operator) and matrix/tensor operators (such as quadrupole and multipole operators).
    using Vectorizer = _Vectorizer;

private:
    // The array that takes care of the storage of the operator's matrix elements.
    StorageArray<QCMatrix<_Scalar>, _Vectorizer> array;


public:
    /*
     *  MARK: Constructors
     */


    /**
     *  Construct a one-electron operator from a(n) (storage) array.
     * 
     *  @param array                A storage array that contains the matrix representations of all the components of this operator.
     */
    SimpleSQOneElectronOperator(const StorageArray<QCMatrix<Scalar>, Vectorizer>& array) :
        array {array} {

        const auto first_dimension = array.elements()[0].dimension();

        for (const auto& parameters : array.elements()) {
            if (parameters.dimension() != first_dimension) {
                throw std::invalid_argument("SimpleSQOneElectronOperator(const StorageArray<QCMatrix<Scalar>, Vectorizer>& array): The dimensions of the matrix representations must be equal.");
            }
        }
    }


    /**
     *  Construct a one-electron operator from one matrix representation.
     * 
     *  @param parameters           The matrix representation of the one-electron integrals/parameters.
     */
    template <typename Z = Vectorizer>
    SimpleSQOneElectronOperator(const QCMatrix<Scalar>& parameters, typename std::enable_if<std::is_same<Z, ScalarVectorizer>::value>::type* = 0) :
        SimpleSQOneElectronOperator(StorageArray<QCMatrix<Scalar>, ScalarVectorizer>({parameters}, ScalarVectorizer())) {}


    /**
     *  Construct a one-electron operator from a set of three matrix representations, for each of the operator's components.
     * 
     *  @param parameters           A set of three matrix representations of the one-electron integrals/parameters, one for each of the operator's components.
     */
    template <typename Z = Vectorizer>
    SimpleSQOneElectronOperator(const std::vector<QCMatrix<Scalar>>& parameters, typename std::enable_if<std::is_same<Z, VectorVectorizer>::value>::type* = 0) :
        SimpleSQOneElectronOperator(StorageArray<QCMatrix<Scalar>, VectorVectorizer>(parameters, VectorVectorizer({3}))) {}


    /**
     *  Construct a one-electron operator with parameters that are zero.
     * 
     *  @param dim          The dimension of the matrix representation of the parameters, i.e. the number of orbitals/sites.
     */
    SimpleSQOneElectronOperator(const size_t dim, const Vectorizer& vectorizer) :
        SimpleSQOneElectronOperator(StorageArray<QCMatrix<Scalar>, Vectorizer> {QCMatrix<Scalar>::Zero(dim, dim), vectorizer}) {}


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
    const std::vector<QCMatrix<Scalar>>& allParameters() const { return this->array.elements(); }

    /**
     *  @return A vector of writable matrix representions of the parameters/integrals of this operator.
     */
    std::vector<QCMatrix<Scalar>>& allParameters() { return this->array.elements(); }

    /**
     *  @param indices      A set of coordinates that accesses this one-electron operator.
     * 
     *  @return A read-only matrix representation of the parameters/integrals of one of the tensor components of this operator.
     */
    template <typename... Indices>
    const QCMatrix<Scalar>& parameters(const Indices&... indices) const { return this->array(indices...); }


    /**
     *  @param indices      A set of coordinates that accesses this one-electron operator.
     * 
     *  @return A writable matrix representation of the parameters/integrals of one of the tensor components of this operator.
     */
    template <typename... Indices>
    QCMatrix<Scalar>& parameters(const Indices&... indices) { return this->array(indices...); }


    // TODO: move to derived class?
    /**
     *  @param indices      A set of coordinates that accesses this one-electron operator.
     * 
     *  @return The component of this one-electron operator that is stored at the given coordinate indices.
     */
    // template <typename... Indices>
    // Element operator()(const Indices&... indices) const {

    //     // Access the underlying array's storage, and wrap the result in a scalar second-quantized operator.
    // }


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
    size_t numberOfOrbitals() const { return this->array.elements()[0].numberOfOrbitals(); /* all the dimensions are the same, this is checked in the constructor */ }


    /*
     *  MARK: Calculations
     */

    /**
     *  @param D                The total 1-DM that represents the wave function
     *
     *  @return The expectation value of this one-electron operator, i.e. the expectation value of all components of the one-electron operator.
     */
    StorageArray<Scalar, Vectorizer> calculateExpectationValue(const OneDM<Scalar>& D) const {

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
};


}  // namespace GQCP
