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
#include "Processing/Properties/RHFElectricalResponseSolver.hpp"

#include "Mathematical/Representation/BlockMatrix.hpp"
#include "QCMethod/HF/RHF.hpp"


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
 *  @param sq_hamiltonian           the Hamiltonian parameters
 * 
 *  @return the parameter response constant (k_p), i.e. the second-order parameter partial derivative of the RHF energy function
 */
SquareMatrix<double> RHFElectricalResponseSolver::calculateParameterResponseConstant(const SQHamiltonian<double>& sq_hamiltonian) const {

    // k_p for RHF is the RHF orbital Hessian
    const auto RHF_orbital_hessian_tensor = calculateRHFOrbitalHessianTensor(sq_hamiltonian, this->N_P);

    return SquareMatrix<double>(RHF_orbital_hessian_tensor.asMatrix());  // give the column-major representation of the tensor as a matrix
}


/**
 *  @param dipole_op                the dipole integrals expressed in an orthonormal orbital basis
 * 
 *  @return the parameter response force (F_p) as an (Nx3)-matrix, i.e. the first-order parameter partial derivative of the perturbation derivative of the RHF energy function
 */
Matrix<double, Dynamic, 3> RHFElectricalResponseSolver::calculateParameterResponseForce(const VectorSQOneElectronOperator<double> dipole_op) const {

    const auto K = dipole_op.dimension();  // number of spatial orbitals

    const auto dim = this->N_P * (K - this->N_P);  // number of non-redundant orbital rotation generators
    Matrix<double, Dynamic, 3> F_p = Matrix<double, Dynamic, 3>::Zero(dim, 3);


    // Loop over the components of the electrical dipole to calculate the response force for every component
    for (size_t m = 0; m < 3; m++) {
        BlockMatrix<double> F_p_m (N_P, K, 0, N_P);  // zero-initialize an object suitable for the representation of virtual-occupied (a,i) quantities

        for (size_t i = 0; i < this->N_P; i++) {
            for (size_t a = this->N_P; a < K; a++) {
                F_p_m(a,i) = 4 * dipole_op.parameters(m)(a,i);  // RHF formula for the parameter (i.e. orbital) response force
            }
        }

        F_p.col(m) = F_p_m.asVector();  // reduce a matrix representation to a column-major vector
    }

    return F_p;
}


}  // namespace GQCP
