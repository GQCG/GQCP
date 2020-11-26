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


#include "Basis/SpinorBasis/RSpinOrbitalBasis.hpp"
#include "Basis/SpinorBasis/USpinOrbitalBasis.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QuantumChemical/Spin.hpp"


namespace GQCP {


/*
 *  MARK: Transformations
 */

/**
 *  Break the recursion of the variadic `transform` method: transform a single object with a basis transformation.
 * 
 *  @tparam TransformationType          The type of the basis transformation.
 *  @tparam TransformableType           The type of the object that should be basis transformed.
 * 
 *  @param T                            The basis transformation.
 *  @param transformable                The object that should be basis transformed.
 */
template <typename TransformationType, typename TransformableType>
void transform(const TransformationType& T, TransformableType& transformable) {
    transformable.transform(T);
}


/**
 *  Basis transform a number of objects with a given transformation.
 * 
 *  @tparam TransformationType          The type of the basis transformation.
 *  @tparam First                       The type of the first object that should be basis transformed.
 *  @tparam Others                      The variadic type of the other objects that should be basis transformed.
 * 
 *  @param T                            The basis transformation.
 *  @param first                        The first object that should be basis transformed.
 *  @param others                       The variadic other objects that should be basis transformed.
 */
template <typename TransformationType, typename First, typename... Others>
void transform(const TransformationType& T, First& first, Others&... others) {

    // Transform the first object, and forward the remaining objects.
    first.transform(T);
    transform(T, others...);
}


/*
 *  MARK: Rotations
 */

/**
 *  Break the recursion of the variadic `rotate` method: rotate a single object with a basis rotation.
 * 
 *  @tparam RotationType                The type of the basis rotation.
 *  @tparam RotatableType               The type of the object that should be basis rotated.
 * 
 *  @param U                            The basis rotation.
 *  @param rotatable                    The object that should be basis rotated.
 */
template <typename RotationType, typename RotatableType>
void rotate(const RotationType& U, RotatableType& rotatable) {
    rotatable.rotate(U);
}


/**
 *  Basis rotate a number of objects with a given transformation.
 * 
 *  @tparam RotationType                The type of the basis rotated.
 *  @tparam First                       The type of the first object that should be basis rotated.
 *  @tparam Others                      The variadic type of the other objects that should be basis rotated.
 * 
 *  @param U                            The basis rotation.
 *  @param first                        The first object that should be basis transformed.
 *  @param others                       The variadic other objects that should be basis transformed.
 */
template <typename RotationType, typename First, typename... Others>
void rotate(const RotationType& U, First& first, Others&... others) {

    // Rotate the first object, and forward the remaining objects.
    first.rotate(U);
    rotate(U, others...);
}


}  // namespace GQCP
