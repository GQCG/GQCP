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


namespace GQCP {


/**
 *  A second-quantized one-electron operator whose parameters are described by a single tensor.
 * 
 *  This class is used as a base class for `RSQTwoElectronOperator` and `GSQTwoElectronOperator`, since they both admit parameter representations using a single tensor, as opposed to `USQTwoElectronOperator`, which uses separate alpha- and beta- tensors. The word 'simple' is used here as an antonym for 'compound'.
 * 
 *  @tparam _Scalar                 The scalar type used for a single parameter/matrix element: real or complex.
 *  @tparam _Vectorizer             The type of the vectorizer that relates a one-dimensional storage of tensors to the tensor structure of two-electron operators. This distinction is carried over from SimpleSQOneElectronOperator.
 *  @tparam _DerivedOperator        The type of the operator that derives from this class, enabling CRTP and compile-time polymorphism.
 */
template <typename _Scalar, typename _Vectorizer, typename _DerivedOperator>
class SimpleSQTwoElectronOperator {};


}  // namespace GQCP
