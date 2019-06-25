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
#include "Properties/RHFElectricalResponseSolver.hpp"

#include "RHF/RHF.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param N_P          the number of electron pairs
 */
RHFElectricalResponseSolver::RHFElectricalResponseSolver(const size_t N_P) :
    N_P (N_P)
{}



/*
 *  PUBLIC OVERRIDDEN METHODS
 */

/**
 *  @param ham_par                  the Hamiltonian parameters
 * 
 *  @return the parameter response constant (k_p), i.e. the second-order parameter partial derivative of the RHF energy function
 */
SquareMatrix<double> RHFElectricalResponseSolver::calculateParameterResponseConstant(const HamiltonianParameters<double>& ham_par) const {

    // k_p for RHF is the RHF orbital Hessian
    const auto RHF_orbital_hessian_tensor = calculateRHFOrbitalHessianTensor(ham_par, this->N_P);

    // We have now calculated the full KxKxKxK tensor (K: number of spatial orbitals)
    // We should reduce the tensor dimension to only include the non-redundant rotations (ignore zero elements) and transform to a Hessian matrix

    return RHF_orbital_hessian_tensor.pairWiseReduce(this->N_P, 0, this->N_P, 0);  // compound both indices (ai, bj) into (m, n), knowing that virtual indices start at N_P
}


/**
 *  @param dipole_integrals         the dipole integrals in an orthonormal orbital basis
 * 
 *  @return the parameter response force (F_p), i.e. the first-order parameter partial derivative of the perturbation derivative of the RHF energy function
 */
Matrix<double, Dynamic, 3> RHFElectricalResponseSolver::calculateParameterResponseForce(const std::array<OneElectronOperator<double>, 3>& dipole_integrals) const {

    const auto K = dipole_integrals[0].get_K();


    const auto dim = this->N_P * (K - this->N_P);  // number of non-redundant orbital rotation generators
    Matrix<double, Dynamic, 3> F_p = Matrix<double, Dynamic, 3>::Zero(dim, 3);


    for (size_t m = 0; m < 3; m++) {  // m loops over the components of the electrical dipole

        // Get only the virtual-occupied dipole integrals
        // FIXME: this code should be gone after the representation refactor
        Matrix<double> mu_m (K, this->N_P);
        mu_m.setZero();
        for (size_t i = 0; i < this->N_P; i++) {
            for (size_t a = this->N_P; a < K; a++) {
                mu_m(a,i) = dipole_integrals[m](a,i);
            }
        }

        // mu_m is a KxK matrix representation of the m-th component of the electrical electron dipole moment
        // Put the virtual-occupied mu_m(a,i) elements into the response force vector
        F_p.col(m) = mu_m.pairWiseReduce(this->N_P, 0);
    }

    return 4 * F_p;
}


}  // namespace GQCP
