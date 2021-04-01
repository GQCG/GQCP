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


#include "Mathematical/Functions/CartesianGTO.hpp"
#include "Mathematical/Functions/EvaluableLinearCombination.hpp"
#include "Mathematical/Functions/Function.hpp"
#include "Mathematical/Representation/Matrix.hpp"
#include "Utilities/aliases.hpp"


namespace GQCP {


/**
 *  One of the elements of the second-quantized field-free charge current density operator.
 *
 *  @param _Scalar              The scalar type of the expansion coefficients of the spatial orbitals that underlie this current density matrix element.
 *  @param _Primitive           The type of primitive that underlies this current density matrix element.
 */
template <typename _Scalar, typename _Primitive>
class CurrentDensityMatrixElement:
    public Function<Vector<complex, 3>, Vector<double, 3>> {
public:
    // The scalar type of the expansion coefficients of the spatial orbitals that underlie this current density matrix element.
    using Scalar = _Scalar;

    // The type of primitive that underlies this current density matrix element.
    using Primitive = _Primitive;

    // The type of basis function that underlies this current density matrix element.
    using BasisFunction = EvaluableLinearCombination<double, Primitive>;

    // The type of spatial orbital that underlies this current density matrix element.
    using SpatialOrbital = EvaluableLinearCombination<Scalar, BasisFunction>;


private:
    SpatialOrbital phi_p;
    SpatialOrbital phi_q;

    SpatialOrbitalGradient dphi_p;
    SpatialOrbitalGradient dphi_q;


public:
    /*
     *  MARK: `Function` behavior
     */

    /**
     *  Evaluate this current density matrix element for the given argument.
     * 
     *  @param r                The point at which the current density matrix element should be evaluated.
     *
     *  @return The function value for the given argument.
     */
    Vector<complex, 3> operator()(const Vector<double, 3>& r) const override {

        return 0.5_ii * (phi_p(r) * dphi_q(r) - dphi_p(r) * phi_q(r));
    }
};


}  // namespace GQCP
