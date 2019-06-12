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
#ifndef GQCP_AP1ROGLAGRANGIANNEWTONORBITALOPTIMIZER_HPP
#define GQCP_AP1ROGLAGRANGIANNEWTONORBITALOPTIMIZER_HPP


#include "OrbitalOptimization/QCMethodNewtonOrbitalOptimizer.hpp"

#include "Geminals/AP1roGGeminalCoefficients.hpp"


namespace GQCP {


class AP1roGLagrangianNewtonOrbitalOptimizer : public QCMethodNewtonOrbitalOptimizer {
private:
    size_t N_P;  // the number of electron pairs
    AP1roGGeminalCoefficients G;  // the current geminal coefficients
    AP1roGVariables multipliers;  // the current Lagrangian multipliers


public:
    // CONSTRUCTORS

    /**
     *  @param N_P              the number of electron pairs
     *  @param G                the initial guess for the AP1roG gemial coefficients
     *  @param oo_options       the options for orbital optimization
     */
    AP1roGLagrangianNewtonOrbitalOptimizer(size_t N_P, const AP1roGGeminalCoefficients& G, const OrbitalOptimizationOptions& oo_options);

    /**
     *  @param N_P              the number of electron pairs
     *  @param K                the number of spatial orbitals
     *  @param oo_options       the options for orbital optimization
     * 
     *  The initial guess for the geminal coefficients is zero
     */
    AP1roGLagrangianNewtonOrbitalOptimizer(size_t N_P, size_t K, const OrbitalOptimizationOptions& oo_options);


    // OVERRIDDEN PUBLIC METHODS

    /**
     *  Prepare this object (i.e. the context for the orbital optimization algorithm) to be able to check for convergence in this Newton-based orbital optimizer
     * 
     *  In the case of this uncoupled AP1roG Lagrangian orbital optimizer, the PSEs are re-solved in every iteration using the current orbitals
     */
    void prepareDMCalculation(const HamiltonianParameters<double>& ham_par) override;

    /**
     *  Prepare this object (i.e. the context for the orbital optimization algorithm) to be able to calculate the new rotation matrix in this Newton-based orbital optimizer
     */
    void prepareQCMethodNewtonSpecificRotationMatrixCalculation(const HamiltonianParameters<double>& ham_par) override {}

    /**
     *  @return the current 1-DM
     */
    OneRDM<double> calculate1RDM() const override;

    /**
     *  @return the current 2-DM
     */
    TwoRDM<double> calculate2RDM() const override;

    /**
     *  Use gradient and Hessian information to determine a new direction for the 'full' orbital rotation generators kappa. Note that a distinction is made between 'free' generators, i.e. those that are calculated from the gradient and Hessian information and the 'full' generators, which also include the redundant parameters (that can be set to zero). The 'full' generators are used to calculate the total rotation matrix using the matrix exponential
     * 
     *  @param ham_par      the current Hamiltonian parameters
     * 
     *  @return the new full set orbital generators, including the redundant parameters
     */
    OrbitalRotationGenerators calculateNewFullOrbitalGenerators(const HamiltonianParameters<double>& ham_par) const override;
};


}  // namespace GQCP


#endif  // GQCP_AP1ROGLAGRANGIANNEWTONORBITALOPTIMIZER_HPP