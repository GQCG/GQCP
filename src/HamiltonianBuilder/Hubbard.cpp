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
#include "HamiltonianBuilder/Hubbard.hpp"


namespace GQCP {

/*
 *  PRIVATE METHODS
 */

/**
 *  Evaluate the one-electron operators for alpha or beta and store the result in a matrix-vector product or a matrix, depending on the method passed
 *
 *  @param fock_space_target        the Fock space that is used as a target, i.e. that is evaluated
 *  @param fock_space_fixed         the Fock space that is not evaluated
 *  @param target_is_major          whether or not the evaluated component is the major index
 *  @param hamiltonian_parameters   the Hubbard Hamiltonian parameters
 *  @param method                   the used method: constructHamiltonian() or matrixVectorProduct()
 */
void Hubbard::oneOperatorModule(FockSpace& fock_space_target, FockSpace& fock_space_fixed, bool target_is_major, const HamiltonianParameters& hamiltonian_parameters, const PassToMethod& method) {

    size_t K = fock_space_target.get_K();
    size_t N = fock_space_target.get_N();
    size_t dim = fock_space_target.get_dimension();
    size_t dim_fixed = fock_space_fixed.get_dimension();

    size_t fixed_intervals;
    size_t target_interval;

    // If the target is major, then the interval for the non-target (or fixed component) is 1
    // while the the target (major) intervals after each fixed (minor) iteration, thus at the dimension of the fixed component.
    if (target_is_major) {
        fixed_intervals = 1;
        target_interval = dim_fixed;

    // vice versa, if the target is not major, its own interval is 1,
    // and the fixed component intervals at the targets dimension.
    } else {
        fixed_intervals = dim;
        target_interval = 1;
    }

    ONV onv = fock_space_target.get_ONV(0);  // onv with address 0
    for (size_t I = 0; I < dim; I++) {  // I loops over all the addresses of the onv
        for (size_t e1 = 0; e1 < N; e1++) {  // e1 (electron 1) loops over the (number of) electrons
            size_t p = onv.get_occupied_index(e1);  // retrieve the index of a given electron

            // remove the weight from the initial address I, because we annihilate
            size_t address = I - fock_space_target.get_vertex_weights(p, e1 + 1);
            // The e2 iteration counts the amount of encountered electrons for the creation operator
            // We only consider greater addresses than the initial one (because of symmetry)
            // Hence we are only required to start counting from the annihilated electron (e1)
            size_t e2 = e1;

            // Starting from e2, the sign is always 1
            // Because the initial sign of the annihilation will also be evaluated for the creation
            // 1*1 or -1*-1 is always 1.
            int sign_e2 = 1;

            // Test whether next orbital is occupied, until we reach unoccupied orbital
            while (e2 < N - 1 && onv.get_occupied_index(e2 + 1) - onv.get_occupied_index(e2) == 1) {
                // Shift the address for the electrons encountered after the annihilation but before the creation
                // Their currents weights are no longer correct, the corresponding weights can be calculated
                // initial weight can be found in the addressing scheme, on the index of the orbital (row) and electron count (column)
                // since e2 starts at the annihilated position, the first shifted electron is at e2's position + 1, (given the while loop condition this is also (e2+1)'s position)
                // The nature of the addressing scheme requires us the add 1 to the electron count (because we start with the 0'th electron
                // And for the initial weight we are at an extra electron (before the annihilation) hence the difference in weight is:
                // the new weight at (e2+1) position (row) and e2+1 (column) - the old weight at  (e2+1) position (row) and e2+2 (column)
                address += fock_space_target.get_vertex_weights(onv.get_occupied_index(e2) + 1, e2 + 1) - fock_space_target.get_vertex_weights(onv.get_occupied_index(e2) + 1, e2 + 2);
                e2++;  // adding occupied orbitals to the electron count
                sign_e2 *= -1;  // skipping over non-annihilated electrons will cause a phase change
            }
            size_t q = onv.get_occupied_index(e2) + 1;
            e2++;
            while (q < K) {
                size_t J = address + fock_space_target.get_vertex_weights(q, e2);

                // address has been calculated, update accordingly and at all instances of the fixed component
                for (size_t I_fixed = 0; I_fixed < dim_fixed; I_fixed++){
                    double val = sign_e2 * hamiltonian_parameters.get_h()(p, q);
                    method(I * target_interval + I_fixed * fixed_intervals, J * target_interval + I_fixed * fixed_intervals, val);
                    method(J * target_interval + I_fixed * fixed_intervals, I * target_interval + I_fixed * fixed_intervals, val);
                }

                // go to the next orbital
                q++;

                // if we encounter an occupied orbital, perform the shift, and test whether the following orbitals are occupied (or not)
                // then proceed to set q to the next non-occupied orbital.
                if (e2 < N && q == onv.get_occupied_index(e2)) {
                    address += fock_space_target.get_vertex_weights(q, e2) - fock_space_target.get_vertex_weights(q, e2 + 1);
                    sign_e2 *= -1;
                    while (e2 < N - 1 && onv.get_occupied_index(e2 + 1) - onv.get_occupied_index(e2) == 1) {
                        // see previous
                        address += fock_space_target.get_vertex_weights(onv.get_occupied_index(e2) + 1, e2 + 1) - fock_space_target.get_vertex_weights(onv.get_occupied_index(e2) + 1, e2 + 2);
                        e2++;
                        sign_e2 *= -1;
                    }
                    q = onv.get_occupied_index(e2) + 1;
                    e2++;
                }
            }  //  (creation)

        } // e1 loop (annihilation)

        // Prevent last permutation
        if (I < dim - 1) {
            fock_space_target.setNext(onv);
        }
    }

}


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
 *  @param hamiltonian_parameters       the Hubbard Hamiltonian parameters in an orthonormal orbital basis
 *
 *  @return the Hubbard Hamiltonian matrix
 */
Eigen::MatrixXd Hubbard::constructHamiltonian(const HamiltonianParameters& hamiltonian_parameters) {
    auto K = hamiltonian_parameters.get_h().get_dim();
    if (K != this->fock_space.get_K()) {
        throw std::invalid_argument("Basis functions of the Fock space and hamiltonian_parameters are incompatible.");
    }

    FockSpace fock_space_alpha = fock_space.get_fock_space_alpha();
    FockSpace fock_space_beta = fock_space.get_fock_space_beta();

    auto dim = fock_space.get_dimension();

    Eigen::MatrixXd result_matrix = Eigen::MatrixXd::Zero(dim, dim);
    result_matrix += this->calculateDiagonal(hamiltonian_parameters).asDiagonal();

    // We pass to a matrix and create the corresponding lambda function
    PassToMethod addToMatrix = [&result_matrix](size_t I, size_t J, double value) { result_matrix(I, J) += value; };

    // perform one electron evaluations, one for the alpha component and one for the beta component.
    // In our case alpha will be major and thus when alpha is the "target" (the operators evaluated)
    // the "target_is_major" will be set to true, when beta is the target it will be set to false
    this->oneOperatorModule(fock_space_alpha, fock_space_beta, true, hamiltonian_parameters, addToMatrix);
    this->oneOperatorModule(fock_space_beta, fock_space_alpha, false, hamiltonian_parameters, addToMatrix);

    return result_matrix;
}


/**
 *  @param hamiltonian_parameters       the Hubbard Hamiltonian parameters in an orthonormal orbital basis
 *  @param x                            the vector upon which the Hubbard Hamiltonian acts
 *  @param diagonal                     the diagonal of the Hubbard Hamiltonian matrix
 *
 *  @return the action of the Hubbard Hamiltonian on the coefficient vector
 */
Eigen::VectorXd Hubbard::matrixVectorProduct(const HamiltonianParameters& hamiltonian_parameters, const Eigen::VectorXd& x, const Eigen::VectorXd& diagonal) {

    auto K = hamiltonian_parameters.get_h().get_dim();
    if (K != this->fock_space.get_K()) {
        throw std::invalid_argument("Basis functions of the Fock space and hamiltonian_parameters are incompatible.");
    }

    FockSpace fock_space_alpha = fock_space.get_fock_space_alpha();
    FockSpace fock_space_beta = fock_space.get_fock_space_beta();

    Eigen::VectorXd matvec = diagonal.cwiseProduct(x);

    // We pass to a the matvec and create the corresponding lambda function
    PassToMethod addToMatvec = [&matvec, &x](size_t I, size_t J, double value) { matvec(I) += value * x(J); };

    // perform one electron evaluations, one for the alpha component and one for the beta component.
    // In our case alpha will be major and thus when alpha is the "target" (the operators evaluated)
    // the "target_is_major" will be set to true, when beta is the target it will be set to false
    this->oneOperatorModule(fock_space_alpha, fock_space_beta, true, hamiltonian_parameters, addToMatvec);
    this->oneOperatorModule(fock_space_beta, fock_space_alpha, false, hamiltonian_parameters, addToMatvec);

    return matvec;
}


/**
 *  @param hamiltonian_parameters       the Hubbard Hamiltonian parameters in an orthonormal orbital basis
 *
 *  @return the diagonal of the matrix representation of the Hubbard Hamiltonian
 */
Eigen::VectorXd Hubbard::calculateDiagonal(const HamiltonianParameters &hamiltonian_parameters) {

    auto K = hamiltonian_parameters.get_h().get_dim();
    if (K != this->fock_space.get_K()) {
        throw std::invalid_argument("Basis functions of the Fock space and hamiltonian_parameters are incompatible.");
    }

    FockSpace fock_space_alpha = fock_space.get_fock_space_alpha();
    FockSpace fock_space_beta = fock_space.get_fock_space_beta();

    auto dim_alpha = fock_space_alpha.get_dimension();
    auto dim_beta = fock_space_beta.get_dimension();
    auto dim = fock_space.get_dimension();

    // Diagonal contributions
    Eigen::VectorXd diagonal = Eigen::VectorXd::Zero(dim);

    ONV onv_alpha = fock_space_alpha.get_ONV(0);
    ONV onv_beta =  fock_space_beta.get_ONV(0);
    for (size_t Ia = 0; Ia < dim_alpha; Ia++) {  // Ia loops over addresses of alpha onvs
        fock_space_beta.set(onv_beta, 0);
        for (size_t Ib = 0; Ib < dim_beta; Ib++) {  // Ib loops over addresses of beta onvs
            size_t address = Ia * dim_beta + Ib;
            // find all double occupations
            Vectoru occupations = onv_alpha.findMatchingOccupations(onv_beta);
            for (size_t p : occupations){
                diagonal(address) += hamiltonian_parameters.get_g()(p,p,p,p);
            }
            if (Ib < dim_beta - 1) {  // prevent last permutation to occur
                fock_space_beta.setNext(onv_beta);
            }
        }  // beta address (Ib) loop
        if (Ia < dim_alpha - 1) {  // prevent last permutation to occur
            fock_space_alpha.setNext(onv_alpha);
        }
    }  // alpha address (Ia) loop
    return diagonal;

}


/*
 *  HELPER METHODS
 */

/**
 *  Generate the upper triagonal as a vector for a Hubbard lattice
 *
 *  @param A        the adjacency matrix that represents the allowed interaction between sites
 *  @param t        the one-electron hopping interaction parameter
 *  @param U        the two-electron interaction parameter
 *
 *  @return the upper triagonal as a vector of the hopping matrix generated from the adjacency matrix and the Hubbard parameters t and U
 */
Eigen::VectorXd generateUpperTriagonal(Eigen::MatrixXd A, double t, double U) {

    // Check if the hopping matrix is represented as a square matrix
    if (A.cols() != A.rows()) {
        throw std::invalid_argument("The adjacency matrix has to be represented as a square matrix.");
    }

    size_t K = A.cols();
    size_t length = (K * (K+1))/2;  // formula for the length of the triagonal (in one vector)
    Eigen::VectorXd triagonal = Eigen::VectorXd::Zero(length);

    size_t index = 0;
    for (size_t i = 0; i < K; i++) {
        for (size_t j = i; j < K; j++) {
            if (i == j) {
                triagonal(index) = U;
            } else {
                triagonal(index) = t * A(i,j);
            }
            index++;
        }
    }

    return triagonal;
}


}  // namespace GQCP
