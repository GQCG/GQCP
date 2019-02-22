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
#include "geminals/BaseAP1roGSolver.hpp"


namespace GQCP {


/*
 * CONSTRUCTORS
 */

/**
 *  @param N_P          the number of electrons
 *  @param ham_par      Hamiltonian parameters in an orthonormal orbital basis
 *  @param G            the initial guess for the AP1roG gemial coefficients
 */
BaseAP1roGSolver::BaseAP1roGSolver(size_t N_P, const HamiltonianParameters<double>& ham_par, const AP1roGGeminalCoefficients& G) :
    K (ham_par.get_K()),
    ham_par (ham_par),
    N_P (N_P),
    geminal_coefficients (G)
{}

/**
 *  @param N_P          the number of electrons
 *  @param ham_par      Hamiltonian parameters in an orthonormal orbital basis
 *
 *  The initial guess for the geminal coefficients is zero
 */
BaseAP1roGSolver::BaseAP1roGSolver(size_t N_P, const HamiltonianParameters<double>& ham_par) :
    BaseAP1roGSolver(N_P, ham_par, AP1roGGeminalCoefficients(N_P, ham_par.get_K()))
{}


/**
 *  @param molecule     the molecule used for the AP1roG calculation
 *  @param ham_par      Hamiltonian parameters in an orthonormal orbital basis
 *  @param G            the initial guess for the AP1roG gemial coefficients
 */
BaseAP1roGSolver::BaseAP1roGSolver(const Molecule& molecule, const HamiltonianParameters<double>& ham_par, const AP1roGGeminalCoefficients& G) :
    BaseAP1roGSolver(molecule.get_N()/2, ham_par, G)
{
    // Check if we have an even number of electrons
    if ((molecule.get_N() % 2) != 0) {
        throw std::invalid_argument("The given number of electrons is odd.");
    }
}


/**
 *  @param molecule     the molecule used for the AP1roG calculation
 *  @param ham_par      Hamiltonian parameters in an orthonormal orbital basis
 *
 *  The initial guess for the geminal coefficients is zero
 */
BaseAP1roGSolver::BaseAP1roGSolver(const Molecule& molecule, const HamiltonianParameters<double>& ham_par) :
    BaseAP1roGSolver(molecule, ham_par, AP1roGGeminalCoefficients(molecule.get_N()/2, ham_par.get_K()))
{}



/*
 *  DESTRUCTOR
 */

BaseAP1roGSolver::~BaseAP1roGSolver() {}


}  // namespace GQCP
