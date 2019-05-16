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
#ifndef ChemicalMatrix_hpp
#define ChemicalMatrix_hpp


#include "math/SquareMatrix.hpp"


namespace GQCP {


/**
 *  An extension of the SquareMatrix class with methods of quantum chemical context
 *
 *  @tparam _Scalar      the scalar type
 */
template<typename _Scalar>
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
     *  Public Methods
     */

    /**
     *  In-place basis transform of the "chemical" matrix
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
    void basisTransform(const SquareMatrix<TransformationScalar> &T) {
        *this = Self(T.adjoint() * (*this) * T);  // this has no aliasing issues (https://eigen.tuxfamily.org/dox/group__TopicAliasing.html)
    }
};


}  // namespace GQCP


#endif /* ChemicalMatrix_hpp */
