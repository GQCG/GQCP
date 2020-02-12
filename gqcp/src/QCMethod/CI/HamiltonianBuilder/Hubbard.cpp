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
 *  @param onv_basis       the full spin-resolved ONV basis
 */
Hubbard::Hubbard(const SpinResolvedONVBasis& onv_basis) :
    onv_basis (onv_basis)
{}



/*
 *  PUBLIC METHODS
 */

/**
 *  @param hubbard_hamiltonian              the Hubbard model Hamiltonian
 *
 *  @return the Hubbard Hamiltonian matrix
 */
SquareMatrix<double> Hubbard::constructHamiltonian(const HubbardHamiltonian<double>& hubbard_hamiltonian) const {

    const auto K = hubbard_hamiltonian.numberOfLatticeSites();
    if (K != this->onv_basis.get_K()) {
        throw std::invalid_argument("Hubbard::constructHamiltonian(const HubbardHamiltonian<double>&): The number of spatial orbitals of this ONV basis and the number of lattice sites for the Hubbard Hamiltonian are incompatible.");
    }


    return this->onv_basis.evaluateOperatorDense(hubbard_hamiltonian.core(), false) + SquareMatrix<double>(this->calculateDiagonal(hubbard_hamiltonian).asDiagonal());
}


/**
 *  @param hubbard_hamiltonian              the Hubbard model Hamiltonian
 *  @param x                                the vector upon which the Hubbard Hamiltonian acts
 *  @param diagonal                         the diagonal of the Hubbard Hamiltonian matrix
 *
 *  @return the action of the Hubbard model Hamiltonian on the coefficient vector
 */
VectorX<double> Hubbard::matrixVectorProduct(const HubbardHamiltonian<double>& hubbard_hamiltonian, const VectorX<double>& x, const VectorX<double>& diagonal) const {

    const auto K = hubbard_hamiltonian.numberOfLatticeSites();
    if (K != this->onv_basis.get_K()) {
        throw std::invalid_argument("Hubbard::constructHamiltonian(const HubbardHamiltonian<double>&): The number of spatial orbitals of this ONV basis and the number of lattice sites for the Hubbard Hamiltonian are incompatible.");
    }


    // Set up ONV bases for alpha and beta
    const SpinUnresolvedONVBasis onv_basis_alpha = onv_basis.get_fock_space_alpha();
    const auto dim_alpha = onv_basis_alpha.get_dimension();

    const SpinUnresolvedONVBasis onv_basis_beta = onv_basis.get_fock_space_beta();
    const auto dim_beta = onv_basis_beta.get_dimension();


    // Calculate the Hubbard matrix-vector product, which is the sum of the alpha- and beta one-electron contributions plus the diagonal (which contains the two-electron contributions).
    VectorX<double> matvec = diagonal.cwiseProduct(x);

    Eigen::Map<Eigen::MatrixXd> matvecmap (matvec.data(), dim_beta, dim_alpha);
    Eigen::Map<const Eigen::MatrixXd> xmap (x.data(), dim_beta, dim_alpha);

    const auto H_alpha = onv_basis_alpha.evaluateOperatorSparse(hubbard_hamiltonian.core(), false);  // false: no diagonal contributions
    const auto H_beta = onv_basis_beta.evaluateOperatorSparse(hubbard_hamiltonian.core(), false);  // false: no diagonal contributions

    matvecmap += xmap * H_alpha + H_beta * xmap;
    return matvec;
}


/**
 *  @param hubbard_hamiltonian              the Hubbard model Hamiltonian
 *
 *  @return the diagonal of the matrix representation of the Hubbard model Hamiltonian
 */
VectorX<double> Hubbard::calculateDiagonal(const HubbardHamiltonian<double>& hubbard_hamiltonian) const {

    const auto K = hubbard_hamiltonian.numberOfLatticeSites();
    if (K != this->onv_basis.get_K()) {
        throw std::invalid_argument("Hubbard::constructHamiltonian(const HubbardHamiltonian<double>&): The number of spatial orbitals of this ONV basis and the number of lattice sites for the Hubbard Hamiltonian are incompatible.");
    }


    // Set up ONV bases for alpha and beta
    const SpinUnresolvedONVBasis onv_basis_alpha = onv_basis.get_fock_space_alpha();
    const auto dim_alpha = onv_basis_alpha.get_dimension();

    SpinUnresolvedONVBasis onv_basis_beta = onv_basis.get_fock_space_beta();
    const auto dim_beta = onv_basis_beta.get_dimension();

    const auto dim = onv_basis.get_dimension();


    // Calculate the diagonal contributions resulting from the two-electron on-site interactions by iterating over all ONVs.
    VectorX<double> diagonal = VectorX<double>::Zero(dim);

    const auto& H = hubbard_hamiltonian.hoppingMatrix();


    SpinUnresolvedONV onv_alpha = onv_basis_alpha.makeONV(0);
    SpinUnresolvedONV onv_beta = onv_basis_beta.makeONV(0);
    for (size_t Ia = 0; Ia < dim_alpha; Ia++) {  // Ia loops over the addresses of alpha ONVs
        onv_basis_beta.transformONV(onv_beta, 0);

        for (size_t Ib = 0; Ib < dim_beta; Ib++) {  // Ib loops over addresses of beta ONVs
            const size_t address = Ia * dim_beta + Ib;  // compound address in the spin-resolved ONV basis

            // There is a contribution for all orbital indices p that are occupied both in the alpha- and beta ONV.
            std::vector<size_t> occupations = onv_alpha.findMatchingOccupations(onv_beta);
            for (const auto& p : occupations) {
                diagonal(address) += H(p,p);  // the two-electron (on-site repulsion) contributions are on the diagonal of the hopping matrix
            }

            if (Ib < dim_beta - 1) {  // prevent the last permutation from occurring
                onv_basis_beta.setNextONV(onv_beta);
            }
        }  // beta address (Ib) loop

        if (Ia < dim_alpha - 1) {  // prevent the last permutation from occurring
            onv_basis_alpha.setNextONV(onv_alpha);
        }
    }  // alpha address (Ia) loop
    return diagonal;
}



}  // namespace GQCP
