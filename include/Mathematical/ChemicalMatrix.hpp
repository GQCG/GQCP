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
#pragma once


#include "Mathematical/SquareMatrix.hpp"
#include "typedefs.hpp"


namespace GQCP {


/**
 *  An extension of the SquareMatrix class with methods of quantum chemical context
 *
 *  @tparam _Scalar      the scalar type
 */
template <typename _Scalar>
class ChemicalMatrix : public SquareMatrix<_Scalar> {
public:
    using Scalar = _Scalar;

    using Base = SquareMatrix<Scalar>;
    using Self = ChemicalMatrix<Scalar>;


public:

    /*
     *  CONSTRUCTORS
     */

    using SquareMatrix<Scalar>::SquareMatrix;  // use base constructors


    /*
     *  GETTERS
     */
    size_t get_K() const { return this->cols(); };


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return the dimension of this matrix representation of the parameters, i.e. the number of orbitals/sites
     */
    size_t dimension() const {
        return this->cols();
    }

    /**
     *  Basis transform this "chemical" matrix
     *
     *  @tparam TransformationScalar        the type of scalar used for the transformation matrix
     *
     *  @param T    the transformation matrix between the old and the new orbital basis, it is used as
     *      b' = b T ,
     *   in which the basis functions are collected as elements of a row vector b
     *
     *  Note that in order to use these transformation formulas, the multiplication between TransformationScalar and Scalar should be 'enabled'. See LinearCombination.hpp for an example
     */
    template <typename TransformationScalar = Scalar>
    auto basisTransform(const SquareMatrix<TransformationScalar> &T) const -> ChemicalMatrix<product_t<Scalar, TransformationScalar>> {

        using ResultScalar = product_t<Scalar, TransformationScalar>;
        return ChemicalMatrix<ResultScalar>(T.adjoint() * (*this) * T);
    }


    /**
     *  In-place rotate the matrix representation of the one-electron operator using a unitary Jacobi rotation matrix constructed from the Jacobi rotation parameters. Note that this function is only available for real (double) matrix representations
     *
     *  @param jacobi_rotation_parameters       the Jacobi rotation parameters (p, q, angle) that are used to specify a Jacobi rotation: we use the (cos, sin, -sin, cos) definition for the Jacobi rotation matrix. See transform() for how the transformation matrix between the two bases should be represented
     */
    template <typename Z = Scalar>
    enable_if_t<std::is_same<Z, double>::value, ChemicalMatrix<Z>> basisRotate(const JacobiRotationParameters& jacobi_rotation_parameters) {

        auto M_copy = *this;

        auto p = jacobi_rotation_parameters.get_p();
        auto q = jacobi_rotation_parameters.get_q();
        auto angle = jacobi_rotation_parameters.get_angle();

        double c = std::cos(angle);
        double s = std::sin(angle);


        // Use Eigen's Jacobi module to apply the Jacobi rotations directly (cfr. T.adjoint() * M * T)
        Eigen::JacobiRotation<double> jacobi (c, s);

        M_copy.applyOnTheLeft(p, q, jacobi.adjoint());
        M_copy.applyOnTheRight(p, q, jacobi);

        return M_copy;
    }
};


}  // namespace GQCP
