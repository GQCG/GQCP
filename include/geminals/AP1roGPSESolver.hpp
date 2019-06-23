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
#ifndef AP1roGPSESolver_hpp
#define AP1roGPSESolver_hpp


#include "Geminals/BaseAP1roGSolver.hpp"


namespace GQCP {


/**
 *  A class that is able to solve the AP1roG PSE equations
 */
class AP1roGPSESolver : public BaseAP1roGSolver {
public:
    // CONSTRUCTORS
    using BaseAP1roGSolver::BaseAP1roGSolver;  // inherit base constructors


    // PUBLIC METHODS

    /**
     *  @param G        the AP1roG geminal coefficients
     *  @param i        the subscript for the coordinate function
     *  @param a        the superscript for the coordinate function
     *  @param j        the subscript for the geminal coefficient
     *  @param b        the superscript for the geminal coefficient
     *
     *  @return the Jacobian element with compound indices (i,a) and (j,b) at the given geminal coefficients
     */
    double calculateJacobianElement(const AP1roGGeminalCoefficients& G, const size_t i, const size_t a, const size_t j, const size_t b) const;

    /**
     *  @param G        the AP1roG geminal coefficients
     *
     *  @return the Jacobian at the given geminal coefficients
     */
    SquareMatrix<double> calculateJacobian(const AP1roGGeminalCoefficients& G) const;

    /**
     *  @param G        the AP1roG geminal coefficients
     *  @param i        the subscript for the coordinate function
     *  @param a        the superscript for the coordinate function
     *
     *  @return the coordinate function with given indices (i,a) at the given geminal coefficients
     */
    double calculateCoordinateFunction(const AP1roGGeminalCoefficients& G, const size_t i, const size_t a) const;

    /**
     *  @param G        the AP1roG geminal coefficients
     *
     *  @return the vector of coordinate functions at the given geminal coefficients
     */
    VectorX<double> calculateCoordinateFunctions(const AP1roGGeminalCoefficients& G) const;

    /**
     *  Set up and solve the projected Schrödinger equations for AP1roG
     */
    void solve() override;
};


}  // namespace GQCP



#endif /* AP1roGPSESolver_hpp */
