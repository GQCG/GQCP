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


    SimpleSQOneElectronOperator(const StorageArray<QCMatrix<_Scalar>, _Vectorizer>& array) :
        array {array} {

        // Check if the array's matrices have equal dimension.
    }


    template <typename Z = Vectorizer>
    SimpleSQOneElectronOperator(const QCMatrix<Scalar>& parameters, typename std::enable_if<std::is_same<Z, ScalarVectorizer>::value>::type* = 0) :
        SimpleSQOneElectronOperator(StorageArray<QCMatrix<Scalar>, ScalarVectorizer>({parameters}, ScalarVectorizer())) {}


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
};


}  // namespace GQCP
