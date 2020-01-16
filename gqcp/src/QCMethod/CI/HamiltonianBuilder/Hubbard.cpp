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
#include "QCMethod/CI/HamiltonianBuilder/Hubbard.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param fock_space       the full alpha and beta product Fock space
 */
Hubbard::Hubbard(const ProductFockSpace& fock_space) :
    HamiltonianBuilder(),
    fock_space(fock_space)
{}



/*
 *  OVERRIDDEN PUBLIC METHODS
 */

/**
 *  @param sq_hamiltonian       the Hubbard Hamiltonian expressed in an orthonormal basis
 *
 *  @return the Hubbard Hamiltonian matrix
 */
SquareMatrix<double> Hubbard::constructHamiltonian(const SQHamiltonian<double>& sq_hamiltonian) const {
    
    auto K = sq_hamiltonian.dimension();
    if (K != this->fock_space.get_K()) {
        throw std::invalid_argument("Hubbard::constructHamiltonian(SQHamiltonian<double>): Basis functions of the Fock space and sq_hamiltonian are incompatible.");
    }

    return this->fock_space.evaluateOperatorDense(sq_hamiltonian.core(), false) + SquareMatrix<double>(this->calculateDiagonal(sq_hamiltonian).asDiagonal());
}


/**
 *  @param sq_hamiltonian               the Hubbard Hamiltonian expressed in an orthonormal basis
 *  @param x                            the vector upon which the Hubbard Hamiltonian acts
 *  @param diagonal                     the diagonal of the Hubbard Hamiltonian matrix
 *
 *  @return the action of the Hubbard Hamiltonian on the coefficient vector
 */
VectorX<double> Hubbard::matrixVectorProduct(const SQHamiltonian<double>& sq_hamiltonian, const VectorX<double>& x, const VectorX<double>& diagonal) const {

    auto K = sq_hamiltonian.dimension();
    if (K != this->fock_space.get_K()) {
        throw std::invalid_argument("Hubbard::matrixVectorProduct(SQHamiltonian<double>, VectorX<double>, VectorX<double>): Basis functions of the Fock space and sq_hamiltonian are incompatible.");
    }

    FockSpace fock_space_alpha = fock_space.get_fock_space_alpha();
    FockSpace fock_space_beta = fock_space.get_fock_space_beta();

    auto dim_alpha = fock_space_alpha.get_dimension();
    auto dim_beta = fock_space_beta.get_dimension();

    VectorX<double> matvec = diagonal.cwiseProduct(x);

    Eigen::Map<Eigen::MatrixXd> matvecmap(matvec.data(), dim_beta, dim_alpha);
    Eigen::Map<const Eigen::MatrixXd> xmap(x.data(), dim_beta, dim_alpha);

    auto beta_hamiltonian = fock_space_beta.evaluateOperatorSparse(sq_hamiltonian.core(), false);
    auto alpha_hamiltonian = fock_space_alpha.evaluateOperatorSparse(sq_hamiltonian.core(), false);

    matvecmap += xmap * alpha_hamiltonian + beta_hamiltonian * xmap;

    return matvec;
}


/**
 *  @param sq_hamiltonian       the Hubbard Hamiltonian expressed in an orthonormal basis
 *
 *  @return the diagonal of the matrix representation of the Hubbard Hamiltonian
 */
VectorX<double> Hubbard::calculateDiagonal(const SQHamiltonian<double>& sq_hamiltonian) const {

    const auto K = sq_hamiltonian.dimension();
    if (K != this->fock_space.get_K()) {
        throw std::invalid_argument("Hubbard::calculateDiagonal(SQHamiltonian<double>): Basis functions of the Fock space and sq_hamiltonian are incompatible.");
    }

    FockSpace fock_space_alpha = fock_space.get_fock_space_alpha();
    FockSpace fock_space_beta = fock_space.get_fock_space_beta();

    const auto dim_alpha = fock_space_alpha.get_dimension();
    const auto dim_beta = fock_space_beta.get_dimension();
    const auto dim = fock_space.get_dimension();
    const auto& g = sq_hamiltonian.twoElectron().parameters();

    // Diagonal contributions
    VectorX<double> diagonal = VectorX<double>::Zero(dim);

    ONV onv_alpha = fock_space_alpha.makeONV(0);
    ONV onv_beta = fock_space_beta.makeONV(0);
    for (size_t Ia = 0; Ia < dim_alpha; Ia++) {  // Ia loops over addresses of alpha onvs
        fock_space_beta.transformONV(onv_beta, 0);

        for (size_t Ib = 0; Ib < dim_beta; Ib++) {  // Ib loops over addresses of beta onvs
            size_t address = Ia * dim_beta + Ib;
            // find all double occupations
            Vectoru occupations = onv_alpha.findMatchingOccupations(onv_beta);
            for (size_t p : occupations) {
                diagonal(address) += g(p,p,p,p);
            }
            if (Ib < dim_beta - 1) {  // prevent last permutation to occur
                fock_space_beta.setNextONV(onv_beta);
            }
        }  // beta address (Ib) loop
        if (Ia < dim_alpha - 1) {  // prevent last permutation to occur
            fock_space_alpha.setNextONV(onv_alpha);
        }
    }  // alpha address (Ia) loop
    return diagonal;

}



}  // namespace GQCP
