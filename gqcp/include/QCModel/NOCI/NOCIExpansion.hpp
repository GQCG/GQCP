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


#include "Basis/NonOrthogonalBasis/GNonOrthogonalStateBasis.hpp"
#include "Basis/NonOrthogonalBasis/RNonOrthogonalStateBasis.hpp"
#include "Basis/NonOrthogonalBasis/UNonOrthogonalStateBasis.hpp"
#include "Mathematical/Representation/Matrix.hpp"
#include "Utilities/aliases.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/dynamic_bitset.hpp>


namespace GQCP {


/**
 *  A class that represents an expansion inside a non-orthogonal basis.
 *
 *  @tparam _Scalar                       The scalar type of the expansion coefficients: real or complex.
 *  @tparam _NonOrthogonalBasis           The type of non-orthogonal basis.
 */
template <typename _Scalar, typename _NonOrthogonalBasis>
class NOCIExpansion {
public:
    // The scalar type of the expansion coefficients: real or complex.
    using Scalar = _Scalar;

    // The type of the non-orthogonal basis.
    using NonOrthogonalBasis = _NonOrthogonalBasis;


private:
    // The non-orthpgonal basis with respect to which the coefficients are defined.
    NonOrthogonalBasis non_orthogonal_basis;

    // The expansion coefficients.
    VectorX<Scalar> expansion_coefficients;

    // The complete state coefficients.
    MatrixX<Scalar> state_coefficients;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  Construct a NOCI expansion inside the given non-orthogonal basis, with corresponding expansion coefficients.
     *
     *  @param non_orthogonal_basis            The non-orthogonal basis with respect to which the coefficients are defined.
     *  @param coefficients                    The expansion coefficients.
     *  @param state                           The total state coefficients.
     */
    LinearExpansion(const ONVBasis& onv_basis, const VectorX<Scalar>& coefficients, const MatrixX<Scalar>& state) :
        non_orthogonal_basis {non_orthogonal_basis},
        expansion_coefficients {coefficients},
        state_coefficeints {state} {}
};


}  // namespace GQCP
