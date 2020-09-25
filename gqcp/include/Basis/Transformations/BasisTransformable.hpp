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
 *  An (abstract) interface for types that may be transformed from one orbital basis to another.
 * 
 *  In general, we adopt the convention outlined in (https://gqcg-res.github.io/knowdes/spinor-transformations.html), where the new orbitals' coefficients can be found in the respective **column** of the related transformation matrix.
 * 
 *  @param T        The type that should be basis-transformable. It is given as a template argument, enabling CRTP.
 *  @param TM       The type of the transformation matrix for which the basis transformation should be defined. // TODO: Rename "TM" to "TransformationMatrix"
 */
template <typename T, typename TM>
class BasisTransformable {
public:
    /*
     *  MARK: Pure virtual methods
     */

    /**
     *  Apply the basis transformation and return the result.
     * 
     *  @param transformation_matrix        The type that encapsulates the basis transformation coefficients.
     * 
     *  @return The basis-transformed object.
     */
    // T transformed(const TM& transformation_matrix) const;
    virtual T transformed(const TM& transformation_matrix) const = 0;


    /*
     *  MARK: Implementations enabled by pure virtual methods
     */

    /**
     *  In-place apply the basis transformation.
     * 
     *  @param transformation_matrix        The type that encapsulates the basis transformation coefficients.
     */
    virtual void transform(const TM& transformation_matrix) {
        static_cast<T&>(*this) = this->transformed(transformation_matrix);
    }


    /**
     *  Apply the basis rotation and return the result.
     * 
     *  @param transformation_matrix        The type that encapsulates the basis rotation coefficients.
     * 
     *  @return The basis-rotated object.
     */
    virtual T rotated(const TM& transformation_matrix) const {

        // Check if the given matrix is actually unitary.
        if (!transformation_matrix.isUnitary(1.0e-12)) {
            throw std::invalid_argument("BasisTransformable::rotated(const TM<Scalar>&): The given transformation matrix is not unitary.");
        }

        return this->transformed(transformation_matrix);
    }


    /**
     *  In-place apply the basis rotation.
     * 
     *  @param transformation_matrix        The type that encapsulates the basis rotation coefficients.
     */
    void rotate(const TM& transformation_matrix) {
        static_cast<T&>(*this) = this->rotated(transformation_matrix);
    }
};


}  // namespace GQCP
