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
#include "HamiltonianBuilder/DOCI.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param fock_space       the full Fock space, identical for alpha and beta
 */
DOCI::DOCI(const FockSpace& fock_space) :
    HamiltonianBuilder(),
    fock_space (fock_space)
{}



/*
 *  OVERRIDDEN PUBLIC METHODS
 */

/**
 *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
 *
 *  @return the DOCI Hamiltonian matrix
 */
Eigen::MatrixXd DOCI::constructHamiltonian(const HamiltonianParameters& hamiltonian_parameters) {
    
    auto K = hamiltonian_parameters.get_h().get_dim();
    if (K != this->fock_space.get_K()) {
        throw std::invalid_argument("Basis functions of the Fock space and hamiltonian_parameters are incompatible.");
    }
    size_t dim = this->fock_space.get_dimension();
    Eigen::VectorXd diagonal = calculateDiagonal(hamiltonian_parameters);
    Eigen::MatrixXd result_matrix = Eigen::MatrixXd::Zero(dim, dim);
    size_t N = this->fock_space.get_N();

    // Create the first spin string. Since in DOCI, alpha == beta, we can just treat them as one
    // And multiply all contributions by 2
    ONV onv = this->fock_space.get_ONV(0);  // spin string with address 0

    for (size_t I = 0; I < dim; I++) {  // I loops over all the addresses of the onv

        result_matrix(I, I) += diagonal(I);

        for (size_t e1 = 0; e1 < N; e1++) {  // e1 (electron 1) loops over the (number of) electrons
            size_t p = onv.get_occupied_index(e1);  // retrieve the index of a given electron

            // remove the weight from the initial address I, because we annihilate
            size_t address = I - this->fock_space.get_vertex_weights(p, e1 + 1);
            // The e2 iteration counts the amount of encountered electrons for the creation operator
            // We only consider greater addresses than the initial one (because of symmetry)
            // Hence we are only required to start counting from the annihilated electron (e1)
            size_t e2 = e1;

            size_t q = p + 1;

            // Test whether next orbital is occupied, until we reach unoccupied orbital
            while (e2 < N - 1 && q == onv.get_occupied_index(e2 + 1)) {
                // Shift the address for the electrons encountered after the annihilation but before the creation
                // Their currents weights are no longer correct, the corresponding weights can be calculated
                // initial weight can be found in the addressing scheme, on the index of the orbital (row) and electron count (column)
                // since e2 starts at the annihilated position, the first shifted electron is at e2's position + 1, (given the while loop condition this is also (e2+1)'s position)
                // The nature of the addressing scheme requires us the add 1 to the electron count (because we start with the 0'th electron
                // And for the initial weight we are at an extra electron (before the annihilation) hence the difference in weight is:
                // the new weight at (e2+1) position (row) and e2+1 (column) - the old weight at  (e2+1) position (row) and e2+2 (column)
                address += this->fock_space.get_vertex_weights(q, e2 + 1) - this->fock_space.get_vertex_weights(q, e2 + 2);
                e2++;  // adding occupied orbitals to the electron count
                q++;
            }

            e2++;
            while (q < K) {
                size_t J = address + this->fock_space.get_vertex_weights(q, e2);

                // address has been calculated, since I is only updated in the outer loop
                // it is more efficient to add it to a locally scoped variable and add it to the vector with one access.

                result_matrix(I, J) += hamiltonian_parameters.get_g()(p, q, p, q);
                result_matrix(J, I) += hamiltonian_parameters.get_g()(p, q, p, q);

                // go to the next orbital
                q++;

                // if we encounter an occupied orbital, perform the shift, and test whether the following orbitals are occupied (or not)
                // then proceed to set q to the next non-occupied orbital.
                if (e2 < N && q == onv.get_occupied_index(e2)) {
                    address += this->fock_space.get_vertex_weights(q, e2) - this->fock_space.get_vertex_weights(q, e2 + 1);
                    q++;
                    while (e2 < N - 1 &&  q == onv.get_occupied_index(e2 + 1)) {
                        // see previous
                        address += this->fock_space.get_vertex_weights(q, e2 + 1) - this->fock_space.get_vertex_weights(q, e2 + 2);
                        e2++;
                        q++;
                    }
                    e2++;
                }
            }  //  (creation)

        } // e1 loop (annihilation)

        // Prevent last permutation
        if (I < dim - 1) {
            this->fock_space.setNext(onv);
        }


    }   // address (I) loop


    return result_matrix;
}


/**
 *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
 *  @param x                            the vector upon which the DOCI Hamiltonian acts
 *  @param diagonal                     the diagonal of the DOCI Hamiltonian matrix
 *
 *  @return the action of the DOCI Hamiltonian on the coefficient vector
 */
Eigen::VectorXd DOCI::matrixVectorProduct(const HamiltonianParameters& hamiltonian_parameters, const Eigen::VectorXd& x, const Eigen::VectorXd& diagonal) {
    auto K = hamiltonian_parameters.get_h().get_dim();
    if (K != this->fock_space.get_K()) {
        throw std::invalid_argument("Basis functions of the Fock space and hamiltonian_parameters are incompatible.");
    }
    size_t dim = this->fock_space.get_dimension();
    // Create the first spin string. Since in DOCI, alpha == beta, we can just treat them as one
    // And multiply all contributions by 2

    ONV onv = this->fock_space.get_ONV(0);  // spin string with address
    size_t N = this->fock_space.get_N();

    // Diagonal contributions
    Eigen::VectorXd matvec = diagonal.cwiseProduct(x);


    for (size_t I = 0; I < dim; I++) {  // I loops over all the addresses of the onv

        double double_I = 0;
        double double_J = x(I);

        for (size_t e1 = 0; e1 < N; e1++) {  // e1 (electron 1) loops over the (number of) electrons
            size_t p = onv.get_occupied_index(e1);  // retrieve the index of a given electron

            // remove the weight from the initial address I, because we annihilate
            size_t address = I - this->fock_space.get_vertex_weights(p, e1 + 1);
            // The e2 iteration counts the amount of encountered electrons for the creation operator
            // We only consider greater addresses than the initial one (because of symmetry)
            // Hence we are only required to start counting from the annihilated electron (e1)
            size_t e2 = e1;

            size_t q = p + 1;

            // Test whether next orbital is occupied, until we reach unoccupied orbital
            while (e2 < N - 1 && q == onv.get_occupied_index(e2 + 1)) {
                // Shift the address for the electrons encountered after the annihilation but before the creation
                // Their currents weights are no longer correct, the corresponding weights can be calculated
                // initial weight can be found in the addressing scheme, on the index of the orbital (row) and electron count (column)
                // since e2 starts at the annihilated position, the first shifted electron is at e2's position + 1, (given the while loop condition this is also (e2+1)'s position)
                // The nature of the addressing scheme requires us the add 1 to the electron count (because we start with the 0'th electron
                // And for the initial weight we are at an extra electron (before the annihilation) hence the difference in weight is:
                // the new weight at (e2+1) position (row) and e2+1 (column) - the old weight at  (e2+1) position (row) and e2+2 (column)
                address += this->fock_space.get_vertex_weights(q, e2 + 1) - this->fock_space.get_vertex_weights(q, e2 + 2);
                e2++;  // adding occupied orbitals to the electron count
                q++;
            }

            e2++;
            while (q < K) {
                size_t J = address + this->fock_space.get_vertex_weights(q, e2);

                // address has been calculated, since I is only updated in the outer loop
                // it is more efficient to add it to a locally scoped variable and add it to the vector with one access.

                double_I += hamiltonian_parameters.get_g()(p, q, p, q) * x(J);
                matvec(J) += hamiltonian_parameters.get_g()(p, q, p, q) * double_J;

                // go to the next orbital
                q++;

                // if we encounter an occupied orbital, perform the shift, and test whether the following orbitals are occupied (or not)
                // then proceed to set q to the next non-occupied orbital.
                if (e2 < N && q == onv.get_occupied_index(e2)) {
                    address += this->fock_space.get_vertex_weights(q, e2) - this->fock_space.get_vertex_weights(q, e2 + 1);
                    q++;
                    while (e2 < N - 1 &&  q == onv.get_occupied_index(e2 + 1)) {
                        // see previous
                        address += this->fock_space.get_vertex_weights(q, e2 + 1) - this->fock_space.get_vertex_weights(q, e2 + 2);
                        e2++;
                        q++;
                    }
                    e2++;
                }
            }  //  (creation)

        } // e1 loop (annihilation)

        // Prevent last permutation
        if (I < dim - 1) {
            this->fock_space.setNext(onv);
        }

        matvec(I) += double_I;

    }   // address (I) loop

    return matvec;
}


/**
 *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
 *
 *  @return the diagonal of the matrix representation of the Hamiltonian given @param hamiltonian_parameters
 */
Eigen::VectorXd DOCI::calculateDiagonal(const HamiltonianParameters& hamiltonian_parameters) {
    size_t dim = this->fock_space.get_dimension();
    Eigen::VectorXd diagonal = Eigen::VectorXd::Zero(dim);
    // Create the first spin string. Since in DOCI, alpha == beta, we can just treat them as one
    // And multiply all contributions by 2
    ONV onv = this->fock_space.get_ONV(0);  // onv with address 0

    for (size_t I = 0; I < dim; I++) {  // I loops over addresses of spin strings
        double double_I = 0;
        for (size_t e1 = 0; e1 < this->fock_space.get_N(); e1++) {  // e1 (electron 1) loops over the (number of) electrons
            size_t p = onv.get_occupied_index(e1);  // retrieve the index of the orbital the electron occupies
            double_I += 2 * hamiltonian_parameters.get_h()(p,p) + hamiltonian_parameters.get_g()(p,p,p,p);
            for (size_t e2 = 0; e2 < e1; e2++) {  // e2 (electron 2) loops over the (number of) electrons
                // Since we are doing a restricted summation q<p (and thus e2<e1), we should multiply by 2 since the summand argument is symmetric.
                size_t q = onv.get_occupied_index(e2);  // retrieve the index of the orbital the electron occupies
                double_I += 2 * (2*hamiltonian_parameters.get_g()(p,p,q,q) - hamiltonian_parameters.get_g()(p,q,q,p));
            }  // q or e2 loop
        } // p or e1 loop

        diagonal(I) += double_I;

        // Skip the last permutation
        if (I < dim-1) {
            this->fock_space.setNext(onv);
        }



    }  // address (I) loop
    return diagonal;
}



}  // namespace GQCP
