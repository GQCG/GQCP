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


#include <stdexcept>


namespace GQCP {


/*
 *  MARK: BasisTransformableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `BasisTransformable`.
 * 
 *  @tparam Type        The type that should conform to `BasisTransformable`.
 */
template <typename Type>
struct BasisTransformableTraits {};


/*
 *  MARK: BasisTransformable
 */

/**
 *  An (abstract) interface for types that may be transformed from one orbital basis to another.
 * 
 *  @param Type         The type that should be basis-transformable. It is given as a template argument, enabling CRTP.
 */
template <typename Type>
class BasisTransformable {
public:
    // The type of the transformation for which the basis transformation should be defined.
    using Transformation = typename BasisTransformableTraits<Type>::Transformation;

public:
    /*
     *  MARK: Pure virtual methods
     */

    /**
     *  Apply the basis transformation and return the result.
     * 
     *  @param T            The basis transformation.
     * 
     *  @return The basis-transformed object.
     */
    virtual Type transformed(const Transformation& T) const = 0;


    /*
     *  MARK: Implementations enabled by pure virtual methods
     */

    /**
     *  In-place apply the basis transformation.
     * 
     *  @param T            The basis transformation.
     */
    virtual void transform(const Transformation& T) {
        static_cast<Type&>(*this) = this->transformed(T);
    }


    /**
     *  Apply the basis rotation and return the result.
     * 
     *  @param U            The basis rotation.
     * 
     *  @return The basis-rotated object.
     */
    virtual Type rotated(const Transformation& U) const {

        // Check if the given transformation is actually unitary.
        if (!U.isUnitary(1.0e-12)) {
            throw std::invalid_argument("BasisTransformable::rotated(const Transformation<Scalar>&): The given transformation is not unitary.");
        }

        return this->transformed(U);
    }


    /**
     *  In-place apply the basis rotation.
     * 
     *  @param U            The basis rotation.
     */
    void rotate(const Transformation& U) {
        static_cast<Type&>(*this) = this->rotated(U);
    }
};


}  // namespace GQCP
