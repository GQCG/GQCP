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
Eigen::MatrixXd DOCI::constructHamiltonian(const HamiltonianParameters& hamiltonian_parameters) const {
    
    auto K = hamiltonian_parameters.get_h().get_dim();
    if (K != this->fock_space.get_K()) {
        throw std::invalid_argument("The number of orbitals for the Fock space and Hamiltonian parameters are incompatible.");
    }
    size_t dim = this->fock_space.get_dimension();
    Eigen::VectorXd diagonal = calculateDiagonal(hamiltonian_parameters);
    Eigen::MatrixXd result_matrix = Eigen::MatrixXd::Zero(dim, dim);
    size_t N = this->fock_space.get_N();


    auto end = this->fock_space.end();  // help the compiler by putting this out the for-loop
    for (auto it = this->fock_space.begin(); it != end; ++it) {

        size_t I = it.currentAddress();
        ONV onv = it.currentONV();

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

                result_matrix(I, J) += hamiltonian_parameters.get_g()(p, q, p, q);
                result_matrix(J, I) += hamiltonian_parameters.get_g()(p, q, p, q);

                q++;  // go to the next orbital

                // perform a shift
                this->fock_space.shiftUntilNextUnoccupiedOrbital<1>(onv, address, q, e2);

            }  // (creation)
        } // e1 loop (annihilation)
    }  // Fock space iteration

    return result_matrix;
}


/**
 *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
 *  @param x                            the vector upon which the DOCI Hamiltonian acts
 *  @param diagonal                     the diagonal of the DOCI Hamiltonian matrix
 *
 *  @return the action of the DOCI Hamiltonian on the coefficient vector
 */
Eigen::VectorXd DOCI::matrixVectorProduct(const HamiltonianParameters& hamiltonian_parameters, const Eigen::VectorXd& x, const Eigen::VectorXd& diagonal) const {

    auto K = hamiltonian_parameters.get_h().get_dim();
    if (K != this->fock_space.get_K()) {
        throw std::invalid_argument("The number of orbitals for the Fock space and Hamiltonian parameters are incompatible.");
    }

    size_t N = this->fock_space.get_N();


    // Start with the diagonal contributions
    Eigen::VectorXd matvec = diagonal.cwiseProduct(x);


    auto end = this->fock_space.end();  // help the compiler by putting this out the for-loop
    for (auto it = this->fock_space.begin(); it != end; ++it) {

        size_t I = it.currentAddress();
        ONV onv = it.currentONV();

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

                double_I += hamiltonian_parameters.get_g()(p, q, p, q) * x(J);
                matvec(J) += hamiltonian_parameters.get_g()(p, q, p, q) * double_J;

                q++;  // go to the next orbital

                // perform a shift
                this->fock_space.shiftUntilNextUnoccupiedOrbital<1>(onv, address, q, e2);

            }  // (creation)
        } // e1 loop (annihilation)

        matvec(I) += double_I;
    }  // Fock space iteration

    return matvec;
}


/**
 *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
 *
 *  @return the diagonal of the matrix representation of the DOCI Hamiltonian
 */
Eigen::VectorXd DOCI::calculateDiagonal(const HamiltonianParameters& hamiltonian_parameters) const {
    size_t dim = this->fock_space.get_dimension();
    Eigen::VectorXd diagonal = Eigen::VectorXd::Zero(dim);

    auto end = this->fock_space.end();  // help the compiler by putting this out the for-loop
    for (auto it = this->fock_space.begin(); it != end; ++it) {

        size_t I = it.currentAddress();
        ONV onv = it.currentONV();

        double double_I = 0;

        for (size_t e1 = 0; e1 < this->fock_space.get_N(); e1++) {  // e1 (electron 1) loops over the (number of) electrons
            size_t p = onv.get_occupation_index(e1);  // retrieve the index of the orbital the electron occupies
            double_I += 2 * hamiltonian_parameters.get_h()(p,p) + hamiltonian_parameters.get_g()(p,p,p,p);
            for (size_t e2 = 0; e2 < e1; e2++) {  // e2 (electron 2) loops over the (number of) electrons
                // Since we are doing a restricted summation q<p (and thus e2<e1), we should multiply by 2 since the summand argument is symmetric.
                size_t q = onv.get_occupation_index(e2);  // retrieve the index of the orbital the electron occupies
                double_I += 2 * (2*hamiltonian_parameters.get_g()(p,p,q,q) - hamiltonian_parameters.get_g()(p,q,q,p));
            }  // q or e2 loop
        } // p or e1 loop

        diagonal(I) += double_I;
    }  // Fock space iterations

    return diagonal;
}



}  // namespace GQCP
