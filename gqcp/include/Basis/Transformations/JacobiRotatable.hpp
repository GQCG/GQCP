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


#include "Basis/Transformations/JacobiRotationParameters.hpp"


namespace GQCP {


/**
 *  An (abstract) interface for types that may be transformed from one orbital basis to another, using a Jacobi rotation.
 * 
 *  In general, we adopt the convention outlined in (https://gqcg-res.github.io/knowdes/spinor-transformations.html), where the new orbitals' coefficients can be found in the respective **column** of the related transformation matrix.
 * 
 *  @param T        The type that should be Jacobi-transformable. It is given as a template argument, enabling CRTP.
 */
template <typename T>
class JacobiRotatable {
public:
    /*
     *  MARK: Pure virtual methods
     */

    /**
     *  Apply the Jacobi rotation and return the result.
     * 
     *  @param jacobi_parameters        The Jacobi rotation parameters.
     * 
     *  @return The jacobi-transformed object.
     */
    virtual T rotated(const JacobiRotationParameters& jacobi_parameters) const = 0;


    /*
     *  MARK: Implementations enabled by pure virtual methods
     */

    /**
     *  In-place apply the Jacobi rotation.
     * 
     *  @param jacobi_parameters        The Jacobi rotation parameters.
     */
    void rotate(const JacobiRotationParameters& jacobi_parameters) {
        static_cast<T&>(*this) = this->rotated(jacobi_parameters);
    }
};


}  // namespace GQCP