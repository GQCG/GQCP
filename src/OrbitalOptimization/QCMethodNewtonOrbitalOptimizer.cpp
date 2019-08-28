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
#include "OrbitalOptimization/QCMethodNewtonOrbitalOptimizer.hpp"


namespace GQCP {


/*
 *  PUBLIC OVERRIDDEN METHODS
 */

/**
 *  Prepare this object (i.e. the context for the orbital optimization algorithm) to be able to check for convergence in this Newton-based orbital optimizer for quantum chemical methods
 * 
 *  In the case of this uncoupled DOCI orbital optimizer, the DOCI eigenvalue problem is re-solved in every iteration using the current orbitals
 */
void QCMethodNewtonOrbitalOptimizer::prepareOrbitalDerivativesCalculation(const HamiltonianParameters<double>& ham_par) {
    
    this->prepareDMCalculation(ham_par);  // this should prepare the calculation of the 1- and 2-DMs

    this->D = this->calculate1RDM();
    this->d = this->calculate2RDM();
}



/**
 *  @param ham_par      the current Hamiltonian parameters
 * 
 *  @return the current orbital gradient as a matrix
 */
SquareMatrix<double> QCMethodNewtonOrbitalOptimizer::calculateGradientMatrix(const HamiltonianParameters<double>& ham_par) const {

    // Calculate the gradient from the Fockian matrix
    const auto F = ham_par.calculateFockianMatrix(this->D, this->d);
    return 2 * (F.parameters() - F.parameters().transpose());
}


/**
 *  @param ham_par      the current Hamiltonian parameters
 * 
 *  @return the current orbital Hessian as a tensor
 */
SquareRankFourTensor<double> QCMethodNewtonOrbitalOptimizer::calculateHessianTensor(const HamiltonianParameters<double>& ham_par) const {

    const auto K = ham_par.get_K();

    // Calculate the Hessian from the super Fockian matrix
    const auto G = ham_par.calculateSuperFockianMatrix(this->D, this->d).parameters();
    SquareRankFourTensor<double> hessian_tensor (K);
    hessian_tensor.setZero();

    for (size_t p = 0; p < K; p++) {
        for (size_t q = 0; q < K; q++) {
            for (size_t r = 0; r < K; r++) {
                for (size_t s = 0; s < K; s++) {
                    hessian_tensor(p,q,r,s) = G(p,q,r,s) - G(p,q,s,r) + G(q,p,s,r) - G(q,p,r,s) + G(r,s,p,q) - G(r,s,q,p) + G(s,r,q,p) - G(s,r,p,q);
                }
            }
        }
    }

    return hessian_tensor;
}


}  // namespace GQCP
