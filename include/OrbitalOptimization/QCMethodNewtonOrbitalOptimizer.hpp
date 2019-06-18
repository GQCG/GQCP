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
#ifndef GQCP_QCMETHODNEWTONORBITALOPTIMIZER_HPP
#define GQCP_QCMETHODNEWTONORBITALOPTIMIZER_HPP


#include "OrbitalOptimization/NewtonOrbitalOptimizer.hpp"


namespace GQCP {


/**
 *  An intermediary base class for orbital optimization of quantum chemical methods: they use the 1- and 2-DM to calculate the gradient and Hessian
 */
class QCMethodNewtonOrbitalOptimizer : public NewtonOrbitalOptimizer {
protected:
    OneRDM<double> D;  // spin-summed 1-RDM
    TwoRDM<double> d;  // spin-summed 2-RDM


public:
    // CONSTRUCTORS
    using NewtonOrbitalOptimizer::NewtonOrbitalOptimizer;  // inherit base constructors


    // DESTRUCTOR
    virtual ~QCMethodNewtonOrbitalOptimizer() = default;


    // GETTERS
    const OneRDM<double>& get_D() const { return this->D; }
    const TwoRDM<double>& get_d() const { return this->d; }


    // PUBLIC PURE VIRTUAL METHODS

    /**
     *  Prepare this object (i.e. the context for the orbital optimization algorithm) to be able to calculate the 1- and 2-DMs
     */
    virtual void prepareDMCalculation(const HamiltonianParameters<double>& ham_par) = 0;

    /**
     *  @return the current 1-DM
     */
    virtual OneRDM<double> calculate1RDM() const = 0;

    /**
     *  @return the current 2-DM
     */
    virtual TwoRDM<double> calculate2RDM() const = 0;


    // PUBLIC OVERRIDDEN METHODS

    /**
     *  Prepare this object (i.e. the context for the orbital optimization algorithm) to be able to check for convergence in this Newton-based orbital optimizer for quantum chemical methods
     * 
     *  In the case of this uncoupled DOCI orbital optimizer, the DOCI eigenvalue problem is re-solved in every iteration using the current orbitals
     */
    void prepareOrbitalDerivativesCalculation(const HamiltonianParameters<double>& ham_par) override;

    /**
     *  @param ham_par      the current Hamiltonian parameters
     * 
     *  @return the current orbital gradient as a matrix
     */
    SquareMatrix<double> calculateGradientMatrix(const HamiltonianParameters<double>& ham_par) const override;

    /**
     *  @param ham_par      the current Hamiltonian parameters
     * 
     *  @return the current orbital Hessian as a tensor
     */
    SquareRankFourTensor<double> calculateHessianTensor(const HamiltonianParameters<double>& ham_par) const override;
};


}  // namespace GQCP



#endif  // GQCP_QCMETHODNEWTONORBITALOPTIMIZER_HPP
