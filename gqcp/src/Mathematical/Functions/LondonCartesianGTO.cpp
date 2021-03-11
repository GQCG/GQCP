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


#include "Mathematical/Functions/LondonCartesianGTO.hpp"

#include "Utilities/literals.hpp"


namespace GQCP {


/*
 *  MARK: Constructors
 */

/**
 *  @param B            The homogeneous magnetic field appearing in the London modification.
 *  @param gto          The base Cartesian GTO.
 */
LondonCartesianGTO::LondonCartesianGTO(const HomogeneousMagneticField& B, const CartesianGTO& gto) :
    B {B},
    gto {gto} {}


/*
 *  MARK: Magnetic field
 */

/**
 *  @return The k-vector of the London plane wave, i.e. the value of the vector potential at the Gaussian center.
 */
Vector<double, 3> LondonCartesianGTO::kVector() const {

    const auto& K = this->gto.center();
    return this->magneticField().vectorPotentialAt(K);
}


/*
 *  MARK: Functional evaluation
 */

/**
 *  Evaluate the London prefactor at a given point in space.
 * 
 *  @param r            The point in space.
 * 
 *  @return The London plane wave phase factor at the given point.
 */
complex LondonCartesianGTO::phaseFactor(const Vector<double, 3>& r) const {

    using namespace GQCP::literals;
    return std::exp(-1.0_ii * this->kVector().dot(r));
}


/**
 *  Evaluate the London (Cartesian) GTO at a given point in space.
 * 
 *  @param r        The point in space.
 *
 *  @return The value of the London (Cartesian) GTO at the given point.
 */
complex LondonCartesianGTO::operator()(const Vector<double, 3>& r) const {

    return this->phaseFactor(r) * this->gto(r);
}


}  // namespace GQCP
