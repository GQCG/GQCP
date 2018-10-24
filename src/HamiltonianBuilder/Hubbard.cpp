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

void Hubbard::oneOperatorModule(FockSpace& fock_space_target, FockSpace& fock_space_fixed, bool target_is_alpha, const HamiltonianParameters& hamiltonian_parameters, const PassToMethod& method) {

    size_t K = fock_space_target.get_K();
    size_t N = fock_space_target.get_N();
    size_t dim = fock_space_target.get_dimension();
    size_t dim_fixed = fock_space_fixed.get_dimension();

    size_t major_dim;
    size_t minor_dim;

    if (target_is_alpha) {
        major_dim = 1;
        minor_dim = dim_fixed;
    } else {
        major_dim = dim;
        minor_dim = 1;
    }

    ONV onv = fock_space_target.get_ONV(0);  // onv with address 0
    for (size_t I = 0; I < dim; I++) {  // Ia loops over all the addresses of the alpha onv
        if (I > 0) {
            fock_space_target.setNext(onv);
        }
        int sign_e1 = -1;  // sign associated with the operation
        for (size_t e1 = 0; e1 < N; e1++) {  // e1 (electron 1) loops over the (number of) electrons
            sign_e1 *= -1;  // sign of the operator a_p_alpha
            size_t p = onv.get_occupied_index(e1);
            size_t address = I - fock_space_target.get_vertex_weights(p, e1 + 1);  // remove the weight from the initial address I
            size_t e2 = e1;  // counting electrons in the orbitals of the coupling spin string |J> with indices < q
            int sign_e2 = sign_e1;

            // search for the first unoccupied index
            while (e2 < N - 1 && onv.get_occupied_index(e2 + 1) - onv.get_occupied_index(e2) == 1) {
                address += fock_space_target.get_vertex_weights(onv.get_occupied_index(e2) + 1, e2 + 1) - fock_space_target.get_vertex_weights(onv.get_occupied_index(e2) + 1, e2 + 2);
                e2++;
                sign_e2 *= -1;
            }
            size_t q = onv.get_occupied_index(e2) + 1;
            e2++;
            while (q < K) {
                size_t J = address + fock_space_target.get_vertex_weights(q, e2);

                for (size_t I_fixed = 0; I_fixed < dim_fixed; I_fixed++){
                    double val = sign_e2 * hamiltonian_parameters.get_h().get(p, q);
                    method(I * minor_dim + I_fixed * major_dim, J * minor_dim + I_fixed * major_dim, val);
                    method(J * minor_dim + I_fixed * major_dim, I * minor_dim + I_fixed * major_dim, val);
                }
                q++;
                if (e2 < N && q == onv.get_occupied_index(e2)) {
                    address += fock_space_target.get_vertex_weights(q, e2) - fock_space_target.get_vertex_weights(q, e2 + 1);

                    while (e2 < N - 1 && onv.get_occupied_index(e2 + 1) - onv.get_occupied_index(e2) == 1) {
                        address += fock_space_target.get_vertex_weights(onv.get_occupied_index(e2) + 1, e2 + 1) - fock_space_target.get_vertex_weights(onv.get_occupied_index(e2) + 1, e2 + 2);
                        e2++;
                        sign_e2 *= -1;
                    }
                    q = onv.get_occupied_index(e2) + 1;
                    e2++;
                }
            }  // e2 loop (creation)

        } // e1 loop (annihilation)
    }

}


/*
 *  CONSTRUCTORS
 */

/**
 *  Constructor given a @param hamiltonian_parameters and @param fock_space
 */

Hubbard::Hubbard(const FockSpaceProduct &fock_space) :
        HamiltonianBuilder(),
        fock_space(fock_space) {}



/*
 *  OVERRIDDEN PUBLIC METHODS
 */

/**
 *  @return the Hamiltonian matrix as an Eigen::MatrixXd given @param hamiltonian_parameters
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

    PassToMethod addToMatrix = [&result_matrix](size_t I, size_t J, double value) { result_matrix(I, J) += value; };
    std::cout<<dim<< " this is dim "<<std::endl;
    this->oneOperatorModule(fock_space_alpha, fock_space_beta, true, hamiltonian_parameters, addToMatrix);
    this->oneOperatorModule(fock_space_beta, fock_space_alpha, false, hamiltonian_parameters, addToMatrix);

    return result_matrix;
}


/**
 *  @return the action of the Hamiltonian (@param hamiltonian_parameters and @param diagonal) on the coefficient vector @param x
 */

Eigen::VectorXd Hubbard::matrixVectorProduct(const HamiltonianParameters& hamiltonian_parameters, const Eigen::VectorXd& x, const Eigen::VectorXd& diagonal) {

    auto K = hamiltonian_parameters.get_h().get_dim();
    if (K != this->fock_space.get_K()) {
        throw std::invalid_argument("Basis functions of the Fock space and hamiltonian_parameters are incompatible.");
    }

    FockSpace fock_space_alpha = fock_space.get_fock_space_alpha();
    FockSpace fock_space_beta = fock_space.get_fock_space_beta();

    Eigen::VectorXd matvec = diagonal.cwiseProduct(x);

    PassToMethod addToMatvec = [&matvec, &x](size_t I, size_t J, double value) { matvec(I) += value * x(J); };

    this->oneOperatorModule(fock_space_alpha, fock_space_beta, true, hamiltonian_parameters, addToMatvec);
    this->oneOperatorModule(fock_space_beta, fock_space_alpha, false, hamiltonian_parameters, addToMatvec);

    return matvec;
}



/**
 *  @return the diagonal of the matrix representation of the Hamiltonian given @param hamiltonian_parameters
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
    for (size_t Ia = 0; Ia < dim_alpha; Ia++) {  // Ia loops over addresses of alpha onvs
        ONV onv_beta = fock_space_beta.get_ONV(0);
        for (size_t Ib = 0; Ib < dim_beta; Ib++) {  // Ib loops over addresses of beta onvs

            size_t address = Ia * dim_beta + Ib;
            Vectoru occupations = onv_alpha.findMatchingOccupations(onv_beta);

            for (size_t p : occupations){

                diagonal(address) += hamiltonian_parameters.get_g().get(p,p,p,p);

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




}  // namespace GQCP
