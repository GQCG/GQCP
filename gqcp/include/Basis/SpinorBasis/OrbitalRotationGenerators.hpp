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


#include "Basis/TransformationMatrix.hpp"
#include "Mathematical/Representation/SquareMatrix.hpp"


namespace GQCP {


/**
 *  A class that represents the orbital rotation generators kappa
 */
class OrbitalRotationGenerators {
private:
    size_t number_of_spatial_orbitals;  // the number of spatial orbitals that can be rotated using these orbital rotation generators
    VectorX<double> kappa_vector;       // the strict upper/lower triangle of the kappa matrix


public:
    // CONSTRUCTORS

    /**
     *  @param  kappa_vector        the orbital rotation generators represented as a vector that corresponds to the strict upper/lower triangle of the kappa matrix
     */
    OrbitalRotationGenerators(const VectorX<double>& kappa_vector);

    /**
     *  @param  kappa_vector        the orbital rotation generators represented as the full antisymmetric matrix kappa
     */
    OrbitalRotationGenerators(const SquareMatrix<double>& kappa_matrix);


    // NAMED CONSTRUCTORS

    /**
     *  Construct orbital rotation generators by adding redundant (i.e. 0) occupied-virtual and virtual-virtual generators to the given occupied-occupied generators
     * 
     *  @param o_o_generators       the occupied-occupied orbital rotation generators
     *  @param K                    the total number of spatial orbitals
     * 
     *  @return 'full' orbital rotation generators from the given occupied-occupied generators
     */
    static OrbitalRotationGenerators FromOccOcc(const OrbitalRotationGenerators& o_o_generators, const size_t K);


    // PUBLIC METHODS

    /**
     *  @return the orbital rotation generators as the strict upper/lower triangle of the kappa matrix
     */
    const VectorX<double>& asVector() const { return this->kappa_vector; }

    /**
     *  @return the antisymmetric orbital rotation generator matrix kappa
     */
    SquareMatrix<double> asMatrix() const;

    /**
     *  @return the unitary matrix that corresponds to these orbital rotation generators, i.e. exp(-kappa)
     */
    TransformationMatrix<double> calculateRotationMatrix() const;

    /*
     *  @return the number of spatial orbitals that can be rotated using these orbital rotation generators
     */
    size_t numberOfSpatialOrbitals() const { return this->number_of_spatial_orbitals; }
};


}  // namespace GQCP
