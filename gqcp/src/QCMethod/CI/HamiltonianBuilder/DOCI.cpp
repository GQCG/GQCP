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
#include "QCMethod/CI/HamiltonianBuilder/DOCI.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param fock_space       the full ONV basis, identical for alpha and beta
 */
DOCI::DOCI(const SpinUnresolvedONVBasis& fock_space) :
    HamiltonianBuilder(),
    fock_space (fock_space)
{}



/*
 *  OVERRIDDEN PUBLIC METHODS
 */

/**
 *  @param sq_hamiltonian           the Hamiltonian expressed in an orthonormal basis
 *
 *  @return the DOCI Hamiltonian matrix
 */
SquareMatrix<double> DOCI::constructHamiltonian(const SQHamiltonian<double>& sq_hamiltonian) const {

    const auto K = sq_hamiltonian.core().get_dim();
    if (K != this->fock_space.get_K()) {
        throw std::invalid_argument("DOCI::constructHamiltonian(SQHamiltonian<double>): The number of orbitals for the ONV basis and Hamiltonian are incompatible.");
    }

    const size_t dim = this->fock_space.get_dimension();
    SquareMatrix<double> result_matrix = SquareMatrix<double>::Zero(dim, dim);
    const size_t N = this->fock_space.get_N();

    const auto& g = sq_hamiltonian.twoElectron().parameters();


    // Create the first spin string. Since in DOCI, alpha == beta, we can just treat them as one and multiply all contributions by 2
    VectorX<double> diagonal = calculateDiagonal(sq_hamiltonian);
    SpinUnresolvedONV onv = this->fock_space.makeONV(0);  // spin string with address 0

    for (size_t I = 0; I < dim; I++) {  // I loops over all the addresses of the onv

        result_matrix(I, I) += diagonal(I);

        for (size_t e1 = 0; e1 < N; e1++) {  // e1 (electron 1) loops over the (number of) electrons
            size_t p = onv.get_occupation_index(e1);  // retrieve the index of a given electron

            // Remove the weight from the initial address I, because we annihilate
            size_t address = I - this->fock_space.get_vertex_weights(p, e1 + 1);
            // The e2 iteration counts the amount of encountered electrons for the creation operator
            // We only consider greater addresses than the initial one (because of symmetry)
            // Hence we only count electron after the annihilated electron (e1)
            size_t e2 = e1 + 1;
            size_t q = p + 1;

            // perform a shift
            this->fock_space.shiftUntilNextUnoccupiedOrbital<1>(onv, address, q, e2);

            while (q < K) {
                size_t J = address + this->fock_space.get_vertex_weights(q, e2);

                result_matrix(I, J) += g(p, q, p, q);
                result_matrix(J, I) += g(p, q, p, q);

                q++;  // go to the next orbital

                // perform a shift
                this->fock_space.shiftUntilNextUnoccupiedOrbital<1>(onv, address, q, e2);

            }  // (creation)

        } // e1 loop (annihilation)

        // Prevent last permutation
        if (I < dim - 1) {
            this->fock_space.setNextONV(onv);
        }


    }  // address (I) loop


    return result_matrix;
}


/**
 *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal basis
 *  @param x                            the vector upon which the DOCI Hamiltonian acts
 *  @param diagonal                     the diagonal of the DOCI Hamiltonian matrix
 *
 *  @return the action of the DOCI Hamiltonian on the coefficient vector
 */
VectorX<double> DOCI::matrixVectorProduct(const SQHamiltonian<double>& sq_hamiltonian, const VectorX<double>& x, const VectorX<double>& diagonal) const {

    const auto K = sq_hamiltonian.core().get_dim();
    if (K != this->fock_space.get_K()) {
        throw std::invalid_argument("DOCI::matrixVectorProduct(SQHamiltonian<double>, VectorX<double>, VectorX<double>): The number of orbitals for the ONV basis and Hamiltonians are incompatible.");
    }

    const size_t N = this->fock_space.get_N();
    const size_t dim = this->fock_space.get_dimension();
    const auto& g = sq_hamiltonian.twoElectron().parameters();


    // Create the first spin string. Since in DOCI, alpha == beta, we can just treat them as one and multiply all contributions by 2
    SpinUnresolvedONV onv = this->fock_space.makeONV(0);  // spin string with address

    // Diagonal contributions
    VectorX<double> matvec = diagonal.cwiseProduct(x);

    for (size_t I = 0; I < dim; I++) {  // I loops over all the addresses of the onv

        // double_I and J reduce vector accessing and writing
        double double_I = 0;
        double double_J = x(I);

        for (size_t e1 = 0; e1 < N; e1++) {  // e1 (electron 1) loops over the (number of) electrons
            size_t p = onv.get_occupation_index(e1);  // retrieve the index of a given electron

            // Remove the weight from the initial address I, because we annihilate
            size_t address = I - this->fock_space.get_vertex_weights(p, e1 + 1);
            // The e2 iteration counts the amount of encountered electrons for the creation operator
            // We only consider greater addresses than the initial one (because of symmetry)
            // Hence we only count electron after the annihilated electron (e1)
            size_t e2 = e1 + 1;
            size_t q = p + 1;

            // perform a shift
            this->fock_space.shiftUntilNextUnoccupiedOrbital<1>(onv, address, q, e2);

            while (q < K) {
                size_t J = address + this->fock_space.get_vertex_weights(q, e2);

                double_I += g(p, q, p, q) * x(J);
                matvec(J) += g(p, q, p, q) * double_J;

                q++;  // go to the next orbital

                // perform a shift
                this->fock_space.shiftUntilNextUnoccupiedOrbital<1>(onv, address, q, e2);

            }  // (creation)

        } // e1 loop (annihilation)

        // Prevent last permutation
        if (I < dim - 1) {
            this->fock_space.setNextONV(onv);
        }

        matvec(I) += double_I;

    }  // address (I) loop

    return matvec;
}


/**
 *  @param sq_hamiltonian               the Hamiltonian parameters in an orthonormal orbital basis
 *
 *  @return the diagonal of the matrix representation of the DOCI Hamiltonian
 */
VectorX<double> DOCI::calculateDiagonal(const SQHamiltonian<double>& sq_hamiltonian) const {

    const auto K = sq_hamiltonian.core().get_dim();
    if (K != this->fock_space.get_K()) {
        throw std::invalid_argument("DOCI::calculateDiagonal(SQHamiltonian<double>): Basis functions of the ONV basis and sq_hamiltonian are incompatible.");
    }

    const size_t dim = this->fock_space.get_dimension();
    const auto& h = sq_hamiltonian.core().parameters();
    const auto& g = sq_hamiltonian.twoElectron().parameters();

    VectorX<double> diagonal = VectorX<double>::Zero(dim);
    // Create the first spin string. Since in DOCI, alpha == beta, we can just treat them as one and multiply all contributions by 2
    SpinUnresolvedONV onv = this->fock_space.makeONV(0);  // onv with address 0

    for (size_t I = 0; I < dim; I++) {  // I loops over addresses of spin strings
        double double_I = 0;

        for (size_t e1 = 0; e1 < this->fock_space.get_N(); e1++) {  // e1 (electron 1) loops over the (number of) electrons
            size_t p = onv.get_occupation_index(e1);  // retrieve the index of the orbital the electron occupies
            double_I += 2 * h(p,p) + g(p,p,p,p);

            for (size_t e2 = 0; e2 < e1; e2++) {  // e2 (electron 2) loops over the (number of) electrons
                // Since we are doing a restricted summation q<p (and thus e2<e1), we should multiply by 2 since the summand argument is symmetric.
                size_t q = onv.get_occupation_index(e2);  // retrieve the index of the orbital the electron occupies
                double_I += 2 * (2*g(p,p,q,q) - g(p,q,q,p));
            }  // q or e2 loop
        } // p or e1 loop

        diagonal(I) += double_I;

        // Skip the last permutation
        if (I < dim-1) {
            this->fock_space.setNextONV(onv);
        }

    }  // address (I) loop
    return diagonal;
}


}  // namespace GQCP
