// This file is part of GQCG-gqcp.
//
// Copyright (C) 2017-2019  the GQCG developers
//
// GQCG-gqcp is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GQCG-gqcp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-gqcp.  If not, see <http://www.gnu.org/licenses/>.
//
#include "Basis/SpinorBasis/OrbitalRotationGenerators.hpp"

#include <unsupported/Eigen/MatrixFunctions>


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param  kappa_vector        the orbital rotation generators represented as a vector that corresponds to the strict upper/lower triangle of the kappa matrix
 */
OrbitalRotationGenerators::OrbitalRotationGenerators(const VectorX<double>& kappa_vector) :
    kappa_vector {kappa_vector},
    number_of_spatial_orbitals {strictTriangularRoot(kappa_vector.size())} {}


/**
 *  @param  kappa_vector        the orbital rotation generators represented as the full antisymmetric matrix kappa
 */
OrbitalRotationGenerators::OrbitalRotationGenerators(const SquareMatrix<double>& kappa_matrix) :
    OrbitalRotationGenerators(kappa_matrix.pairWiseStrictReduce()) {}


/*
 *  NAMED CONSTRUCTORS
 */

/**
 *  Construct orbital rotation generators by adding redundant (i.e. 0) occupied-virtual and virtual-virtual generators to the given occupied-occupied generators
 * 
 *  @param o_o_generators       the occupied-occupied orbital rotation generators
 *  @param K                    the number of spatial orbitals
 * 
 *  @return 'full' orbital rotation generators from the given occupied-occupied generators
 */
OrbitalRotationGenerators OrbitalRotationGenerators::FromOccOcc(const OrbitalRotationGenerators& o_o_generators, const size_t K) {

    SquareMatrix<double> kappa_full_matrix = SquareMatrix<double>::Zero(K, K);

    kappa_full_matrix.topLeftCorner(o_o_generators.numberOfSpatialOrbitals(), o_o_generators.numberOfSpatialOrbitals()) = o_o_generators.asMatrix();
    return OrbitalRotationGenerators(kappa_full_matrix);
}


/*
 *  PUBLIC METHODS
 */

/**
 *  @return the antisymmetric orbital rotation generator matrix kappa
 */
SquareMatrix<double> OrbitalRotationGenerators::asMatrix() const {

    const auto kappa_matrix = GQCP::SquareMatrix<double>::FromStrictTriangle(this->kappa_vector);  // lower triangle only
    const GQCP::SquareMatrix<double> kappa_matrix_transpose = kappa_matrix.transpose();
    return kappa_matrix - kappa_matrix_transpose;  // add the antisymmetric component
}


/**
 *  @return the unitary matrix that corresponds to these orbital rotation generators, i.e. exp(-kappa)
 */
TransformationMatrix<double> OrbitalRotationGenerators::calculateRotationMatrix() const {
    return (-this->asMatrix()).exp();
}


}  // namespace GQCP
