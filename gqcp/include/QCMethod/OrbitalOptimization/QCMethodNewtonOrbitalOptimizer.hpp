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


#include "DensityMatrix/Orbital1DM.hpp"
#include "DensityMatrix/Orbital2DM.hpp"
#include "QCMethod/OrbitalOptimization/NewtonOrbitalOptimizer.hpp"


namespace GQCP {


/**
 *  An intermediary base class for orbital optimization of quantum chemical methods: they use the 1- and 2-DM to calculate the gradient and Hessian
 */
class QCMethodNewtonOrbitalOptimizer:
    public NewtonOrbitalOptimizer {
protected:
    Orbital1DM<double> D;  // spin-summed 1-DM
    Orbital2DM<double> d;  // spin-summed 2-DM


public:
    // CONSTRUCTORS
    using NewtonOrbitalOptimizer::NewtonOrbitalOptimizer;  // inherit base constructors


    // DESTRUCTOR

    /**
     *  The default destructor.
     */
    virtual ~QCMethodNewtonOrbitalOptimizer() = default;


    // PUBLIC PURE VIRTUAL METHODS

    /**
     *  @return the current 1-DM
     */
    virtual Orbital1DM<double> calculate1DM() const = 0;

    /**
     *  @return the current 2-DM
     */
    virtual Orbital2DM<double> calculate2DM() const = 0;

    /**
     *  Prepare this object (i.e. the context for the orbital optimization algorithm) to be able to calculate the 1- and 2-DMs
     */
    virtual void prepareDMCalculation(const SQHamiltonian<double>& sq_hamiltonian) = 0;


    // PUBLIC OVERRIDDEN METHODS

    /**
     *  @param sq_hamiltonian      the current Hamiltonian
     * 
     *  @return the current orbital gradient as a matrix
     */
    SquareMatrix<double> calculateGradientMatrix(const SQHamiltonian<double>& sq_hamiltonian) const override;

    /**
     *  @param sq_hamiltonian      the current Hamiltonian
     * 
     *  @return the current orbital Hessian as a tensor
     */
    SquareRankFourTensor<double> calculateHessianTensor(const SQHamiltonian<double>& sq_hamiltonian) const override;

    /**
     *  Prepare this object (i.e. the context for the orbital optimization algorithm) to be able to check for convergence in this Newton-based orbital optimizer for quantum chemical methods.
     */
    void prepareOrbitalDerivativesCalculation(const SQHamiltonian<double>& sq_hamiltonian) override;


    // PUBLIC METHODS

    /**
     *  @return the 1-DM calculated by this orbital optimizer
     */
    const Orbital1DM<double>& oneDM() const { return this->D; }

    /**
     *  @return the 2-DM calculated by this orbital optimizer
     */
    const Orbital2DM<double>& twoDM() const { return this->d; }
};


}  // namespace GQCP
