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


#include "Mathematical/Representation/DenseVectorizer.hpp"


namespace GQCP {


/*
 *  MARK: BaseFQOperator implementation
 */

/**
 *  A base class used to represent first-quantized one-electron operators.
 * 
 *  @tparam _N                  The number of electrons related to this operator: usually 1 or 2.
 *  @tparam _Scalar             The scalar representation of the operator: real or complex.
 *  @tparam _Vectorizer         The type of the vectorizer that relates to the tensor structure of one-electron operators. This allows for a distinction between scalar operators (such as the kinetic energy operator), vector operators (such as the spin operator) and matrix/tensor operators (such as quadrupole and multipole operators).
 */
template <size_t N, typename _Scalar, typename _Vectorizer>
class BaseFQOperator {
public:
    // The number of electrons related to this operator: usually 1 or 2.
    static constexpr auto NumberOfElectrons = N;

    // The scalar representation of the operator: real or complex.
    using Scalar = _Scalar;

    // The type of the vectorizer that relates to the tensor structure of one-electron operators. This allows for a distinction between scalar operators (such as the kinetic energy operator), vector operators (such as the spin operator) and matrix/tensor operators (such as quadrupole and multipole operators).
    using Vectorizer = _Vectorizer;


public:
    /*
     *  MARK: Destructor
     */

    virtual ~BaseFQOperator() {}
};


/*
 *  MARK: BaseScalarFQOperator
 */

/**
 *  A base class for scalar operators. Since the underlying `Vectorizer` is `ScalarVectorizer`, we're sure -at compile time- that the number of components is 1.
 * 
 *  @tparam _N                  The number of electrons related to this operator: usually 1 or 2.
 *  @tparam _Scalar             The scalar representation of the operator: real or complex.
 */
template <size_t N, typename _Scalar>
class BaseScalarFQOperator:
    public BaseFQOperator<N, _Scalar, ScalarVectorizer> {
public:
    /*
     *  MARK: Vectorizer
     */

    // The number of components of the operator.
    static constexpr size_t NumberOfComponents = 1;

    // The vectorizer related to this operator.
    static const ScalarVectorizer vectorizer;
};

// Instantiate the static const vectorizer.
template <size_t N, typename _Scalar>
const ScalarVectorizer BaseScalarFQOperator<N, _Scalar>::vectorizer {};


/*
 *  MARK: Convenience aliases for one-electron operators
 */

// A scalar-like one-electron operator.
template <typename Scalar>
using BaseScalarFQOneElectronOperator = BaseScalarFQOperator<1, Scalar>;

// A vector-like one-electron operator.
template <typename Scalar>
using BaseVectorFQOneElectronOperator = BaseFQOperator<1, Scalar, VectorVectorizer>;

// A matrix-like one-electron operator.
template <typename Scalar>
using BaseMatrixFQOneElectronOperator = BaseFQOperator<1, Scalar, MatrixVectorizer>;

// A tensor-like one-electron operator.
template <typename Scalar, size_t N>
using BaseTensorFQOneElectronOperator = BaseFQOperator<1, Scalar, TensorVectorizer<N>>;


/*
 *  MARK: Convenience aliases for two-electron operators
 */

// A scalar-like two-electron operator.
template <typename Scalar>
using BaseScalarFQTwoElectronOperator = BaseScalarFQOperator<2, Scalar>;

// A vector-like two-electron operator.
template <typename Scalar>
using BaseVectorFQTwoElectronOperator = BaseFQOperator<2, Scalar, VectorVectorizer>;

// A matrix-like two-electron operator.
template <typename Scalar>
using BaseMatrixFQTwoElectronOperator = BaseFQOperator<2, Scalar, MatrixVectorizer>;

// A tensor-like two-electron operator.
template <typename Scalar, size_t N>
using BaseTensorFQTwoElectronOperator = BaseFQOperator<2, Scalar, TensorVectorizer<N>>;


}  // namespace GQCP
