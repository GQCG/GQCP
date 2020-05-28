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


#include "Mathematical/Representation/ImplicitMatrixSlice.hpp"
#include "QCMethod/OrbitalOptimization/QCMethodNewtonOrbitalOptimizer.hpp"
#include "QCModel/Geminals/AP1roGGeminalCoefficients.hpp"


namespace GQCP {


class AP1roGLagrangianNewtonOrbitalOptimizer: public QCMethodNewtonOrbitalOptimizer {
private:
    size_t N_P;  // the number of electron pairs

    double pse_convergence_threshold;         // the threshold used to check for convergence on the geminal coefficients
    size_t pse_maximum_number_of_iterations;  // maximum number of Newton steps that may be used to achieve convergence of the PSEs

    double E;                                   // the electronic energy
    AP1roGGeminalCoefficients G;                // the current geminal coefficients
    ImplicitMatrixSlice<double> m_multipliers;  // the current Lagrangian m_multipliers


public:
    // CONSTRUCTORS

    /**
     *  @param N_P                                      the number of electron pairs
     *  @param K                                        the number of spatial orbitals
     *  @param hessian_modifier                         the modifier functor that should be used when an indefinite Hessian is encountered
     *  @param oo_convergence_threshold                 the threshold used to check for convergence on the orbital gradient
     *  @param oo_maximum_number_of_iterations          the maximum number of orbital rotation iterations that may be used to achieve convergence
     *  @param pse_convergence_threshold                the threshold used to check for convergence on the geminal coefficients
     *  @param pse_maximum_number_of_iterations         the maximum number of Newton steps that may be used to achieve convergence of the PSEs
     *
     *  The initial guess for the geminal coefficients is zero
     */
    AP1roGLagrangianNewtonOrbitalOptimizer(const size_t N_P, const size_t K, std::shared_ptr<BaseHessianModifier> hessian_modifier, const double oo_convergence_threshold = 1.0e-08, const size_t oo_maximum_number_of_iterations = 128, const double pse_convergence_threshold = 1.0e-08, const size_t pse_maximum_number_of_iterations = 128);

    /**
     *  @param G                                        the initial geminal coefficients
     *  @param hessian_modifier                         the modifier functor that should be used when an indefinite Hessian is encountered
     *  @param oo_convergence_threshold                 the threshold used to check for convergence
     *  @param oo_maximum_number_of_iterations          the maximum number of iterations that may be used to achieve convergence
     *  @param pse_convergence_threshold                the threshold used to check for convergence on the geminal coefficients
     *  @param pse_maximum_number_of_iterations         the maximum number of Newton steps that may be used to achieve convergence of the PSEs
     */
    AP1roGLagrangianNewtonOrbitalOptimizer(const AP1roGGeminalCoefficients& G, std::shared_ptr<BaseHessianModifier> hessian_modifier, const double oo_convergence_threshold = 1.0e-08, const size_t oo_maximum_number_of_iterations = 128, const double pse_convergence_threshold = 1.0e-08, const size_t pse_maximum_number_of_iterations = 128);


    // PUBLIC OVERRIDDEN METHODS

    /**
     *  Prepare this object (i.e. the context for the orbital optimization algorithm) to be able to check for convergence in this Newton-based orbital optimizer
     * 
     *  In the case of this uncoupled AP1roG Lagrangian orbital optimizer, the PSEs are re-solved in every iteration using the current orbitals
     */
    void prepareDMCalculation(const SQHamiltonian<double>& sq_hamiltonian) override;

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
     *  @param sq_hamiltonian      the current Hamiltonian
     * 
     *  @return the new full set orbital generators, including the redundant parameters
     */
    OrbitalRotationGenerators calculateNewFullOrbitalGenerators(const SQHamiltonian<double>& sq_hamiltonian) const override;


    // PUBLIC METHODS

    /**
     *  @return the electronic energy calculated by this orbital optimizer
     */
    double electronicEnergy() const { return this->E; }

    /**
     *  @return the geminal coefficients calculated by this orbital optimizer
     */
    const AP1roGGeminalCoefficients& geminalCoefficients() const { return this->G; }

    /**
     *  @return the Lagrange multipliers calculated by this orbital optimizer
     */
    const ImplicitMatrixSlice<double>& multipliers() const { return this->m_multipliers; }
};


}  // namespace GQCP
