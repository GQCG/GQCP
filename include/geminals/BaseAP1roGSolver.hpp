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
#ifndef BaseAP1roGSolver_hpp
#define BaseAP1roGSolver_hpp


#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "Molecule.hpp"
#include "geminals/AP1roGGeminalCoefficients.hpp"


namespace GQCP {


/**
 *  A base class for solvers using the AP1roG wave function
 */
class BaseAP1roGSolver {
protected:
    size_t K;  // the number of spatial orbitals
    size_t N_P;  // the number of electron pairs
    double electronic_energy;  // the converged electronic energy

    AP1roGGeminalCoefficients geminal_coefficients;  // the converged geminal coefficients

    HamiltonianParameters<double> ham_par;


public:
    // CONSTRUCTORS
    /**
     *  @param N_P          the number of electrons
     *  @param ham_par      Hamiltonian parameters in an orthonormal orbital basis
     *  @param G            the initial guess for the AP1roG gemial coefficients
     */
    BaseAP1roGSolver(size_t N_P, const HamiltonianParameters<double>& ham_par, const AP1roGGeminalCoefficients& G);

    /**
     *  @param N_P          the number of electrons
     *  @param ham_par      Hamiltonian parameters in an orthonormal orbital basis
     *
     *  The initial guess for the geminal coefficients is zero
     */
    BaseAP1roGSolver(size_t N_P, const HamiltonianParameters<double>& ham_par);

    /**
     *  @param molecule     the molecule used for the AP1roG calculation
     *  @param ham_par      Hamiltonian parameters in an orthonormal orbital basis
     *  @param G            the initial guess for the AP1roG gemial coefficients
     */
    BaseAP1roGSolver(const Molecule& molecule, const HamiltonianParameters<double>& ham_par, const AP1roGGeminalCoefficients& G);

    /**
     *  @param molecule     the molecule used for the AP1roG calculation
     *  @param ham_par      Hamiltonian parameters in an orthonormal orbital basis
     *
     *  The initial guess for the geminal coefficients is zero
     */
    BaseAP1roGSolver(const Molecule& molecule, const HamiltonianParameters<double>& ham_par);


    // DESTRUCTOR
    virtual ~BaseAP1roGSolver();


    // GETTERS
    double get_electronic_energy() const { return this->electronic_energy; }
    const AP1roGGeminalCoefficients& get_geminal_coefficients() const { return this->geminal_coefficients; }
    const HamiltonianParameters<double>& get_ham_par() const { return this->ham_par; }


    // PUBLIC METHODS
    /**
     *  The actual 'solving' step
     */
    virtual void solve() = 0;
};


}  // namespace GQCP


#endif /* BaseAP1roGSolver_hpp */
