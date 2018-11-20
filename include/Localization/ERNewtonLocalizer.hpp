// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2018  the GQCG developers
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
#ifndef ERNewtonLocalizer_hpp
#define ERNewtonLocalizer_hpp


#include "Localization/BaseERLocalizer.hpp"


namespace GQCP {



/**
 *  A class that localizes a set of orthonormal orbitals according to the maximization of the Edmiston-Ruedenberg localization index. A maximum is found using subsequent Newton steps.
 */
class ERNewtonLocalizer : public BaseERLocalizer {
private:
    // PRIVATE METHODS
    /**
     *  @param ham_par      the Hamiltonian parameters containing the two-electron integrals
     *  @param i            the row of the gradient 'matrix'
     *  @param j            the column of the gradient 'matrix'
     *
     *  @return the element (i,j) of the Edmiston-Ruedenberg localization index gradient
     */
    double calculateGradientElement(const GQCP::HamiltonianParameters& ham_par, size_t i, size_t j) const;

    /**
     *  @param ham_par      the Hamiltonian parameters containing the two-electron integrals
     *
     *  @return the gradient of the Edmiston-Ruedenberg localization index as a matrix
     */
    Eigen::MatrixXd calculateGradient(const GQCP::HamiltonianParameters& ham_par) const;

    /**
     *  @param ham_par      the Hamiltonian parameters containing the two-electron integrals
     *  @param i            the first index of the Hessian 'tensor'
     *  @param j            the second index of the Hessian 'tensor'
     *  @param k            the third index of the Hessian 'tensor'
     *  @param l            the fourth index of the Hessian 'tensor'
     *
     *  @return the element (i,j,k,l) of the Edmiston-Ruedenberg localization index Hessian
     */
    double calculateHessianElement(const GQCP::HamiltonianParameters& ham_par, size_t i, size_t j, size_t k, size_t l) const;

    /**
     *  @param ham_par      the Hamiltonian parameters containing the two-electron integrals
     *
     *  @return the Hessian of the Edmiston-Ruedenberg localization index as a tensor
     */
    Eigen::Tensor<double, 4> calculateHessian(const GQCP::HamiltonianParameters& ham_par) const;


public:
    // CONSTRUCTORS
    /**
     *  @param N_P                              the number of electron pairs
     *  @param threshold                        the threshold for maximization on subsequent localization indices
     *  @param maximum_number_of_iterations     the maximum number of iterations for the localization algorithm
     */
    ERNewtonLocalizer(size_t N_P, double threshold=1.0e-08, size_t maximum_number_of_iterations=128);


    // PUBLIC METHODS
    /**
     *  Localize the Hamiltonian parameters by maximizing the Edmiston-Ruedenberg localization index, using a Newton-based algorithm
     *
     *  @param ham_par      the Hamiltonian parameters that should be localized
     */
    void localize(GQCP::HamiltonianParameters& ham_par) override;
};



}  // namespace GQCP


#endif /* ERNewtonLocalizer_hpp */
