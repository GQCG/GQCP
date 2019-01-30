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


#include "geminals/BaseAP1roGSolver.hpp"


namespace GQCP {


/**
 *  A class that is able to solve the AP1roG PSE equations
 */
class AP1roGPSESolver : public BaseAP1roGSolver {
public:
    // CONSTRUCTORS
    /**
     *  @param N_P          the number of electrons
     *  @param ham_par      Hamiltonian parameters in an orthonormal orbital basis
     *  @param G            the initial guess for the AP1roG gemial coefficients
     */
    AP1roGPSESolver(size_t N_P, const HamiltonianParameters& ham_par, const AP1roGGeminalCoefficients& G);

    /**
     *  @param N_P          the number of electrons
     *  @param ham_par      Hamiltonian parameters in an orthonormal orbital basis
     *
     *  The initial guess for the geminal coefficients is zero
     */
    AP1roGPSESolver(size_t N_P, const HamiltonianParameters& ham_par);

    /**
     *  @param molecule     the molecule used for the AP1roG calculation
     *  @param ham_par      Hamiltonian parameters in an orthonormal orbital basis
     *  @param G            the initial guess for the AP1roG gemial coefficients
     */
    AP1roGPSESolver(const Molecule& molecule, const HamiltonianParameters& ham_par, const AP1roGGeminalCoefficients& G);

    /**
     *  @param molecule     the molecule used for the AP1roG calculation
     *  @param ham_par      Hamiltonian parameters in an orthonormal orbital basis
     *
     *  The initial guess for the geminal coefficients is zero
     */
    AP1roGPSESolver(const Molecule& molecule, const HamiltonianParameters& ham_par);


    // PUBLIC METHODS
    /**
     *  @param G        the AP1roG geminal coefficients
     *  @param i        the subscript for the coordinate function
     *  @param a        the superscript for the coordinate function
     *  @param k        the subscript for the geminal coefficient
     *  @param c        the superscript for the geminal coefficient
     *
     *  @return the Jacobian element with compound indices (i,a) and (k,c) at the given geminal coefficients
     */
    double calculateJacobianElement(const AP1roGGeminalCoefficients& G, size_t i, size_t a, size_t k, size_t c) const;

    /**
     *  @param g        the AP1roG geminal coefficients in row-major vector form
     *
     *  @return the Jacobian at the given geminal coefficients
     */
    Eigen::MatrixXd calculateJacobian(const Eigen::VectorXd& g) const;

    /**
     *  @param G        the AP1roG geminal coefficients
     *  @param i        the subscript for the coordinate function
     *  @param a        the superscript for the coordinate function
     *
     *  @return the coordinate function with given indices (i,a) at the given geminal coefficients
     */
    double calculateCoordinateFunction(const AP1roGGeminalCoefficients& G, size_t i, size_t a) const;

    /**
     *  @param g        the AP1roG geminal coefficients in row-major vector form
     *
     *  @return the vector of coordinate functions at the given geminal coefficients
     */
    Eigen::VectorXd calculateCoordinateFunctions(const Eigen::VectorXd& g) const;

    /**
     *  Set up and solve the projected Schrödinger equations for AP1roG
     */
    void solve() override;
};


}  // namespace GQCP



#endif /* AP1roGPSESolver_hpp */
