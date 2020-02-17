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
 *  @param onv_basis       the full seniority-zero ONV basis
 */
DOCI::DOCI(const SeniorityZeroONVBasis& onv_basis) :
    onv_basis (onv_basis)
{}



/*
 *  PUBLIC METHODS
 */

/**
 *  @param sq_hamiltonian           the Hamiltonian expressed in an orthonormal basis
 *
 *  @return the DOCI Hamiltonian matrix
 */
SquareMatrix<double> DOCI::constructHamiltonian(const SQHamiltonian<double>& sq_hamiltonian) const {

    const auto K = sq_hamiltonian.dimension();  // the number of spatial orbitals

    if (K != this->onv_basis.numberOfSpatialOrbitals()) {
        throw std::invalid_argument("DOCI::constructHamiltonian(const SQHamiltonian<double>&): The number of spatial orbitals for the ONV basis and Hamiltonian are incompatible.");
    }


    // Prepare some variables to be used in the algorithm.
    const size_t N = this->onv_basis.numberOfElectrons();
    const size_t dim = this->onv_basis.dimension();

    SquareMatrix<double> result_matrix = SquareMatrix<double>::Zero(dim, dim);
    const auto diagonal = this->calculateDiagonal(sq_hamiltonian);
    const auto& g = sq_hamiltonian.twoElectron().parameters();


    // Create the first doubly-occupied ONV basis. Since in DOCI, alpha == beta, we can use the proxy ONV basis to treat them as one and multiply all contributions by 2.
    const proxy_onv_basis = this->onv_basis.proxy();

    auto onv = proxy_onv_basis.makeONV(0);  // ONV with address 0
    for (size_t I = 0; I < dim; I++) {  // I loops over all the addresses of the onv

        result_matrix(I, I) += diagonal(I);

        for (size_t e1 = 0; e1 < N; e1++) {  // e1 (electron 1) loops over the number of electrons
            const size_t p = onv.get_occupation_index(e1);  // retrieve the index of the orbital that the electron occupies

            // Remove the weight from the initial address I, because we annihilate
            const size_t address = I - proxy_onv_basis.get_vertex_weights(p, e1 + 1);

            // The e2 iteration counts the number of encountered electrons for the creation operator.
            // We only consider greater addresses than the initial one (because of symmetry), hence we only count electron after the annihilated electron (e1).
            const size_t e2 = e1 + 1;
            size_t q = p + 1;

            // perform a shift  TODO: clarify what shift
            proxy_onv_basis.shiftUntilNextUnoccupiedOrbital<1>(onv, address, q, e2);

            while (q < K) {
                const size_t J = address + proxy_onv_basis.get_vertex_weights(q, e2);

                result_matrix(I, J) += g(p,q,p,q);
                result_matrix(J, I) += g(p,q,p,q);

                q++;  // go to the next orbital

                // perform a shift  TODO: clarify what shift
                proxy_onv_basis.shiftUntilNextUnoccupiedOrbital<1>(onv, address, q, e2);
            }  // creation
        } // e1 loop (annihilation)

        if (I < dim - 1) {  // prevent the last permutation from occurring
            proxy_onv_basis.setNextONV(onv);
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

    const auto K = sq_hamiltonian.dimension();
    if (K != this->onv_basis.numberOfSpatialOrbitals()) {
        throw std::invalid_argument("DOCI::matrixVectorProduct(const SQHamiltonian<double>&, const VectorX<double>&, const VectorX<double>&): The number of spatial orbitals for the ONV basis and Hamiltonian are incompatible.");
    }

    // Prepare some variables to be used in the algorithm.
    const size_t N = this->onv_basis.numberOfElectrons();
    const size_t dim = this->onv_basis.dimension();

    const auto& g = sq_hamiltonian.twoElectron().parameters();


    // Initialize the resulting matrix-vector product from the diagonal contributions.
    VectorX<double> matvec = diagonal.cwiseProduct(x);

    // Create the first doubly-occupied ONV basis. Since in DOCI, alpha == beta, we can use the proxy ONV basis to treat them as one and multiply all contributions by 2.
    const proxy_onv_basis = this->onv_basis.proxy();
    auto onv = proxy_onv_basis.makeONV(0);  // ONV with address 0
    for (size_t I = 0; I < dim; I++) {  // I loops over all the addresses of the onv

        // Using container values of type double reduce the number of times a vector has to be read from/written to
        double value = 0;
        const double x_I = x(I);

        for (size_t e1 = 0; e1 < N; e1++) {  // e1 (electron 1) loops over the (number of) electrons
            const size_t p = onv.get_occupation_index(e1);  // retrieve the index of a given electron

            // Remove the weight from the initial address I, because we annihilate
            size_t address = I - proxy_onv_basis.get_vertex_weights(p, e1 + 1);

            // The e2 iteration counts the number of encountered electrons for the creation operator.
            // We only consider greater addresses than the initial one (because of symmetry), hence we only count electron after the annihilated electron (e1).
            const size_t e2 = e1 + 1;
            size_t q = p + 1;

            // perform a shift  TODO: clarify what shift
            proxy_onv_basis.shiftUntilNextUnoccupiedOrbital<1>(onv, address, q, e2);

            while (q < K) {
                const size_t J = address + proxy_onv_basis.get_vertex_weights(q, e2);

                value += g(p, q, p, q) * x(J);
                matvec(J) += g(p, q, p, q) * x_I;

                q++;  // go to the next orbital

                // perform a shift  TODO: clarify what shift
                proxy_onv_basis.shiftUntilNextUnoccupiedOrbital<1>(onv, address, q, e2);
            }  // creation
        } // e1 loop (annihilation)

        if (I < dim - 1) {  // prevent last permutation
            proxy_onv_basis.setNextONV(onv);
        }

        matvec(I) += value;
    }  // address (I) loop

    return matvec;
}


/**
 *  @param sq_hamiltonian               the Hamiltonian parameters in an orthonormal orbital basis
 *
 *  @return the diagonal of the matrix representation of the DOCI Hamiltonian
 */
VectorX<double> DOCI::calculateDiagonal(const SQHamiltonian<double>& sq_hamiltonian) const {

    const auto K = sq_hamiltonian.dimension();  // number of spatial orbitals

    if (K != this->onv_basis.numberOfSpatialOrbitals()) {
        throw std::invalid_argument("DOCI::calculateDiagonal(const SQHamiltonian<double>&): The number of spatial orbitals for the ONV basis and Hamiltonian are incompatible.");
    }

    // Prepare some variables to be used in the algorithm.
    const size_t dim = this->onv_basis.dimension();
    const auto& h = sq_hamiltonian.core().parameters();
    const auto& g = sq_hamiltonian.twoElectron().parameters();

    VectorX<double> diagonal = VectorX<double>::Zero(dim);


    // Create the first doubly-occupied ONV basis. Since in DOCI, alpha == beta, we can use the proxy ONV basis to treat them as one and multiply all contributions by 2.
    const proxy_onv_basis = this->onv_basis.proxy();
    SpinUnresolvedONV onv = proxy_onv_basis.makeONV(0);  // ONV with address 0

    for (size_t I = 0; I < dim; I++) {  // I loops over addresses of spin strings

        double value = 0;  // to be added to the diagonal

        for (size_t e1 = 0; e1 < proxy_onv_basis.numberOfElectrons(); e1++) {  // e1 (electron 1) loops over the number of electrons
            const size_t p = onv.get_occupation_index(e1);  // retrieve the index of the orbital that the electron occupies
            value += 2 * h(p,p) + g(p,p,p,p);

            for (size_t e2 = 0; e2 < e1; e2++) {  // e2 (electron 2) loops over the number of electrons

                // Since we are doing a restricted summation (e2 < e1), we should multiply by 2 since the summand argument is symmetric.
                const size_t q = onv.get_occupation_index(e2);  // retrieve the index of the orbital the electron occupies
                value += 2 * (2*g(p,p,q,q) - g(p,q,q,p));
            }  // q or e2 loop
        } // p or e1 loop

        diagonal(I) += value;

        if (I < dim-1) {  // prevent the last permutation from occurring
            proxy_onv_basis.setNextONV(onv);
        }
    }  // address (I) loop

    return diagonal;
}


}  // namespace GQCP
