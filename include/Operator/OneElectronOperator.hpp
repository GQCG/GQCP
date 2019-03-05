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
#ifndef GQCP_ONEELECTRONOPERATOR_HPP
#define GQCP_ONEELECTRONOPERATOR_HPP


#include "JacobiRotationParameters.hpp"
#include "math/SquareMatrix.hpp"
#include "Operator.hpp"
#include "math/ScalarFunction.hpp"


namespace GQCP {


/**
 *  A class that represents a one-electron operator in an orbital basis
 *
 *  @tparam _Scalar      the scalar type
 */
template <typename _Scalar>
class OneElectronOperator : public SquareMatrix<_Scalar>, public Operator<OneElectronOperator<_Scalar>> {
public:

    using Scalar = _Scalar;

    using BaseRepresentation = SquareMatrix<Scalar>;
    using Self = OneElectronOperator<Scalar>;


public:

    /*
     *  CONSTRUCTORS
     */

    using SquareMatrix<Scalar>::SquareMatrix;  // use base constructors


    /*
     *  PUBLIC METHODS
     */

    /**
     *  In-place transform the matrix representation of the one-electron operator
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
    void transform(const SquareMatrix<TransformationScalar>& T) {
        *this = OneElectronOperator<Scalar>(T.adjoint() * (*this) * T);  // this has no aliasing issues (https://eigen.tuxfamily.org/dox/group__TopicAliasing.html)
    }


    using Operator<OneElectronOperator<Scalar>>::rotate;  // bring over rotate from the base class


    /**
     *  In-place rotate the matrix representation of the one-electron operator using a unitary Jacobi rotation matrix constructed from the Jacobi rotation parameters. Note that this function is only available for real (double) matrix representations
     *
     *  @param jacobi_rotation_parameters       the Jacobi rotation parameters (p, q, angle) that are used to specify a Jacobi rotation: we use the (cos, sin, -sin, cos) definition for the Jacobi rotation matrix. See transform() for how the transformation matrix between the two bases should be represented
     */
    template<typename Z = Scalar>
    enable_if_t<std::is_same<Z, double>::value> rotate(const JacobiRotationParameters& jacobi_rotation_parameters) {

        auto p = jacobi_rotation_parameters.get_p();
        auto q = jacobi_rotation_parameters.get_q();
        auto angle = jacobi_rotation_parameters.get_angle();

        double c = std::cos(angle);
        double s = std::sin(angle);


        // Use Eigen's Jacobi module to apply the Jacobi rotations directly (cfr. T.adjoint() * M * T)
        Eigen::JacobiRotation<double> jacobi (c, s);

        this->applyOnTheLeft(p, q, jacobi.adjoint());
        this->applyOnTheRight(p, q, jacobi);
    }


    /**
     *  @param x        the vector/point at which the scalar functions should be evaluated
     *
     *  @return a one-electron operator corresponding to the evaluated scalar functions
     *
     *  Note that this function is only available for OneElectronOperators whose Scalar is a derived class of ScalarFunction
     */
    template <typename Z = Scalar>
    enable_if_t<std::is_base_of<ScalarFunction<typename Z::Valued, typename Z::Scalar, Z::Cols>, Z>::value,
    OneElectronOperator<typename Z::Valued>> evaluate(const Vector<typename Z::Scalar, Z::Cols>& x) const {

        Eigen::Matrix<typename Z::Valued, SquareMatrix<Z>::Rows, SquareMatrix<Z>::Cols> result (this->rows(), this->cols());
        auto result_op = SquareMatrix<typename Z::Valued>(result);

        for (size_t i = 0; i < this->rows(); i++) {
            for (size_t j = 0; j < this->cols(); j++) {
                result_op(i,j) = (*this)(i,j).operator()(x);
            }
        }

        return OneElectronOperator<typename Z::Valued>(result_op);
    }
};



}  // namespace GQCP



#endif  // GQCP_ONEELECTRONOPERATOR_HPP
