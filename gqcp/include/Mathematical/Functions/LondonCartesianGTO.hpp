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
#include "Mathematical/Functions/Function.hpp"
#include "Physical/HomogeneousMagneticField.hpp"
#include "Utilities/aliases.hpp"


namespace GQCP {


/**
 *  A London-modified Cartesian GTO.
 *
 *  @note Mathematically speaking, a Cartesian GTO is a complex-valued scalar function taking an Euclidean vector (3D-vector) as argument, which is why we inherit from `Function`.
 *  @note Calling operator() returns the value of the unnormalized London GTO.
 *  @note Contracted London GTOs can be expressed as linear combinations of GTOs: `EvaluableLinearCombination<LondonCartesianGTO>`.
 */
class LondonCartesianGTO:
    public Function<complex, double, 3> {
public:
    // The return type of the call operator, i.e. the valuedness of the scalar function.
    using Valued = complex;


private:
    // The homogeneous magnetic field appearing in the London modification.
    HomogeneousMagneticField B;

    // The base Cartesian GTO.
    CartesianGTO gto;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  @param B            The homogeneous magnetic field appearing in the London modification.
     *  @param gto          The base Cartesian GTO.
     */
    LondonCartesianGTO(const HomogeneousMagneticField& B, const CartesianGTO& gto);


    /*
     *  MARK: Magnetic field
     */

    /**
     *  @return The homogeneous magnetic field appearing in the London modification.
     */
    const HomogeneousMagneticField& magneticField() const { return this->B; }

    /**
     *  @return The k-vector of the London plane wave, i.e. the value of the vector potential at the Gaussian center.
     */
    Vector<double, 3> kVector() const;


    /*
     *  MARK: Function evaluation
     */

    /**
     *  @return The base Cartesian GTO.
     */
    const CartesianGTO& cartesianGTO() const { return this->gto; }

    /**
     *  Evaluate the London prefactor at a given point in space.
     * 
     *  @param r            The point in space.
     * 
     *  @return The London plane wave phase factor at the given point.
     */
    complex phaseFactor(const Vector<double, 3>& r) const;

    /**
     *  Evaluate the London (Cartesian) GTO at a given point in space.
     * 
     *  @param r        The point in space.
     *
     *  @return The value of the London (Cartesian) GTO at the given point.
     */
    complex operator()(const Vector<double, 3>& r) const;
};


}  // namespace GQCP
