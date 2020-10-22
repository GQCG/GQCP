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

#include "Processing/Properties/RHFElectricalResponseSolver.hpp"

#include "Basis/SpinorBasis/OrbitalSpace.hpp"
#include "Mathematical/Representation/ImplicitMatrixSlice.hpp"
#include "QCMethod/HF/RHF/RHF.hpp"
#include "QCModel/HF/RHF.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param N_P          the number of electron pairs
 */
RHFElectricalResponseSolver::RHFElectricalResponseSolver(const size_t N_P) :
    N_P {N_P} {}


/*
 *  PUBLIC OVERRIDDEN METHODS
 */

/**
 *  @param sq_hamiltonian           the Hamiltonian parameters
 * 
 *  @return the parameter response constant (k_p), i.e. the second-order parameter partial derivative of the RHF energy function
 */
SquareMatrix<double> RHFElectricalResponseSolver::calculateParameterResponseConstant(const RSQHamiltonian<double>& sq_hamiltonian) const {

    // k_p for RHF is the RHF orbital Hessian
    const auto RHF_orbital_hessian_tensor = QCModel::RHF<double>::calculateOrbitalHessianTensor(sq_hamiltonian, this->N_P);

    return SquareMatrix<double>(RHF_orbital_hessian_tensor.asMatrix());  // give the column-major representation of the tensor as a matrix
}


/**
 *  @param dipole_op                the dipole integrals expressed in an orthonormal orbital basis
 * 
 *  @return the parameter response force (F_p) as an (Nx3)-matrix, i.e. the first-order parameter partial derivative of the perturbation derivative of the RHF energy function
 */
Matrix<double, Dynamic, 3> RHFElectricalResponseSolver::calculateParameterResponseForce(const VectorRSQOneElectronOperator<double>& dipole_op) const {

    // Create an occupied-virtual orbital space.
    const auto K = dipole_op.numberOfOrbitals();  // number of spatial orbitals

    const auto orbital_space = OrbitalSpace::Implicit({{OccupationType::k_occupied, N_P}, {OccupationType::k_virtual, K - N_P}});  // N_P occupied (spatial) orbitals, K-N_P virtual (spatial) orbitals


    // Zero-initialize the vector representation of the parameter response force F_p
    const auto dim = orbital_space.numberOfExcitations(OccupationType::k_occupied, OccupationType::k_virtual);  // number of non-redundant orbital rotation generators
    Matrix<double, Dynamic, 3> F_p = Matrix<double, Dynamic, 3>::Zero(dim, 3);

    // Loop over the components of the electrical dipole to calculate the response force for every component
    for (size_t m = 0; m < 3; m++) {
        auto F_p_m = orbital_space.initializeRepresentableObjectFor<double>(OccupationType::k_virtual, OccupationType::k_occupied);  // zero-initialize an object suitable for the representation of virtual-occupied (a,i) quantities

        for (const auto& i : orbital_space.indices(OccupationType::k_occupied)) {
            for (const auto& a : orbital_space.indices(OccupationType::k_virtual)) {
                F_p_m(a, i) = 4 * dipole_op.parameters(m)(a, i);  // RHF formula for the parameter (i.e. orbital) response force
            }
        }

        F_p.col(m) = F_p_m.asVector();  // reduce a matrix representation to a column-major vector
    }

    return F_p;
}


}  // namespace GQCP
