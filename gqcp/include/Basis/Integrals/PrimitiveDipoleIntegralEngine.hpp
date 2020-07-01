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
#include "Operator/FirstQuantized/ElectronicDipoleOperator.hpp"


namespace GQCP {


/**
 *  A class that can calculate electronic dipole integrals over primitive Cartesian GTOs.
 */
class PrimitiveDipoleIntegralEngine {
public:
    static constexpr auto Components = ElectronicDipoleOperator::Components;
    using IntegralScalar = ElectronicDipoleOperator::Scalar;


private:
    ElectronicDipoleOperator dipole_operator;  // the dipole operator over which this engine should calculate integrals

    CartesianDirection component;  // the current component of the dipole operator this engine can calculate integrals over


public:
    // CONSTRUCTORS
    /**
     *  Construct a PrimitiveDipoleIntegralEngine from its members.
     * 
     *  @param dipole_operator              the dipole operator over which this engine should calculate integrals
     *  @param component                    the current component of the dipole operator this engine can calculate integrals over
     */
    PrimitiveDipoleIntegralEngine(const ElectronicDipoleOperator& dipole_operator, const CartesianDirection component = CartesianDirection::x);


    // PUBLIC METHODS

    /**
     *  @param left             the left Cartesian GTO (primitive)
     *  @param right            the right Cartesian GTO (primitive)
     * 
     *  @return the dipole integral over the two given primitives
     */
    IntegralScalar calculate(const CartesianGTO& left, const CartesianGTO& right);

    /**
     *  @param alpha            the Gaussian exponent of the left 1-D primitive
     *  @param K                the (directional coordinate of the) center of the left 1-D primitive
     *  @param i                the Cartesian exponent of the left 1-D primitive
     *  @param beta             the Gaussian exponent of the right 1-D primitive
     *  @param L                the (directional coordinate of the) center of the right 1-D primitive
     *  @param j                the Cartesian exponent of the right 1-D primitive
     * 
     *  @return the dipole integral over the two given 1-D primitives
     */
    IntegralScalar calculate1D(const double alpha, const double K, const int i, const double beta, const double L, const int j);

    /**
     *  Prepare this engine's internal state such that it is able to calculate integrals over the given component of the operator.
     * 
     *  @param component                the index of the component of the operator
     * 
     *  @note For the dipole integrals, the components are the x, y, and z-components of the (relative) position operator.
     */
    void prepareStateForComponent(const size_t component);
};


}  // namespace GQCP
