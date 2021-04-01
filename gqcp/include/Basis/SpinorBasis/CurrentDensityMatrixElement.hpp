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
#include "Utilities/complex.hpp"


namespace GQCP {


/**
 *  One of the elements of (on of the components of) the second-quantized field-free charge current density operator.
 *
 *  @param _Scalar              The scalar type of the expansion coefficients of the spatial orbitals that underlie this current density matrix element.
 *  @param _Primitive           The type of primitive that underlies this current density matrix element.
 */
template <typename _Scalar, typename _Primitive>
class CurrentDensityMatrixElement:
    public Function<complex, Vector<double, 3>> {
public:
    // The scalar type of the expansion coefficients of a spatial orbital.
    using Scalar = _Scalar;

    // The type of primitive that underlies this current density matrix element.
    using Primitive = _Primitive;

    // The type of basis function that underlies this current density matrix element.
    using BasisFunction = EvaluableLinearCombination<double, Primitive>;

    // The type of spatial orbital that underlies this current density matrix element.
    using SpatialOrbital = EvaluableLinearCombination<Scalar, BasisFunction>;

    // The type of the derivative of a primitive. The derivative of a Cartesian GTO is a linear combination of Cartesian GTOs.
    using PrimitiveDerivative = EvaluableLinearCombination<double, Primitive>;

    // The type of the derivative of a basis function.
    using BasisFunctionDerivative = EvaluableLinearCombination<double, PrimitiveDerivative>;

    // The type of the derivative of a spatial orbital.
    using SpatialOrbitalDerivative = EvaluableLinearCombination<Scalar, BasisFunctionDerivative>;


private:
    // The first spatial orbital.
    SpatialOrbital phi_p;

    // The second spatial orbital.
    SpatialOrbital phi_q;

    // One of the components of the first spatial orbital gradient.
    SpatialOrbitalDerivative dphi_p;

    // One of the components of the second spatial orbital gradient.
    SpatialOrbitalDerivative dphi_q;


public:
    /*
     *  MARK: Constructors
     */

    /*
     *  The default constructor.
     */
    CurrentDensityMatrixElement() = default;


    /**
     *  @param phi_p            The first spatial orbital.
     *  @param phi_q            The second spatial orbital.
     *  @param dphi_p           One of the components of the first spatial orbital gradient.
     *  @param dphi_q           One of the components of the second spatial orbital gradient.
     */
    CurrentDensityMatrixElement(const SpatialOrbital& phi_p, const SpatialOrbital& phi_q, const SpatialOrbitalDerivative& dphi_p, const SpatialOrbitalDerivative& dphi_q) :
        phi_p {phi_p},
        phi_q {phi_q},
        dphi_p {dphi_p},
        dphi_q {dphi_q} {}


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
    complex operator()(const Vector<double, 3>& r) const override {

        using namespace GQCP::literals;

        return 0.5_ii * (conj(phi_p(r)) * dphi_q(r) -
                         conj(dphi_p(r)) * phi_q(r));
    }
};


}  // namespace GQCP
