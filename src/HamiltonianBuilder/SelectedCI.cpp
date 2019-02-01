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
#include "HamiltonianBuilder/SelectedCI.hpp"


namespace GQCP {

/*
 *  PRIVATE METHODS
 */

/**
 *  Evaluate the hamiltonian elements
 *
 *  @param hamiltonian_parameters   the orthogonal Hamiltonian parameters
 *  @param method                   the used method: constructHamiltonian() or matrixVectorProduct()
 */
void SelectedCI::evaluateHamiltonian(const HamiltonianParameters& hamiltonian_parameters, const PassToMethod& method) const {

    size_t dim = fock_space.get_dimension();
    size_t N_a = fock_space.get_N_alpha();
    size_t N_b = fock_space.get_N_beta();

    for (size_t I = 0; I < dim; I++) {  // loop over all addresses (1)
        Configuration configuration_I = this->fock_space.get_configuration(I);
        ONV alpha_I = configuration_I.onv_alpha;
        ONV beta_I = configuration_I.onv_beta;

        // Calculate the off-diagonal elements, by going over all other ONVs
        for (size_t J = I+1; J < dim; J++) {

            Configuration configuration_J = this->fock_space.get_configuration(J);
            ONV alpha_J = configuration_J.onv_alpha;
            ONV beta_J = configuration_J.onv_beta;

            if ((alpha_I.countNumberOfDifferences(alpha_J) == 2) && (beta_I.countNumberOfDifferences(beta_J) == 0)) {

                // Find the orbitals that are occupied in one string, and aren't in the other
                size_t p = alpha_I.findDifferentOccupations(alpha_J)[0];  // we're sure that there is only 1 element in the std::vector<size_t>
                size_t q = alpha_J.findDifferentOccupations(alpha_I)[0];  // we're sure that there is only 1 element in the std::vector<size_t>

                // Calculate the total sign
                int sign = alpha_I.operatorPhaseFactor(p) * alpha_J.operatorPhaseFactor(q);

                double value = hamiltonian_parameters.get_h()(p,q);

                method(I, J, sign*value);
                method(J, I, sign*value);

                for (size_t r = 0; r < K; r++) {  // r loops over spatial orbitals

                    if (alpha_I.isOccupied(r) && alpha_J.isOccupied(r)) {  // r must be occupied on the left and on the right
                        if ((p != r) && (q != r)) {  // can't create or annihilate the same orbital

                            double value = hamiltonian_parameters.get_g()(p,q,r,r)
                                         - hamiltonian_parameters.get_g()(r,q,p,r)
                                         - hamiltonian_parameters.get_g()(p,r,r,q)
                                         + hamiltonian_parameters.get_g()(r,r,p,q);

                            method(I, J, sign*value);
                            method(J, I, sign*value);
                        }
                    }

                    if (beta_I.isOccupied(r)) {  // beta_I == beta_J from the previous if-branch

                        double value = hamiltonian_parameters.get_g()(p,q,r,r)
                                       + hamiltonian_parameters.get_g()(q,p,r,r);

                        method(I, J, sign*value);
                        method(J, I, sign*value);
                    }
                }
            }

            // 0 electron excitations in alpha, 1 in beta
            if ((alpha_I.countNumberOfDifferences(alpha_J) == 0) && (beta_I.countNumberOfDifferences(beta_J) == 2)) {


                // Find the orbitals that are occupied in one string, and aren't in the other
                size_t p = beta_I.findDifferentOccupations(beta_J)[0];  // we're sure that there is only 1 element in the std::vector<size_t>
                size_t q = beta_J.findDifferentOccupations(beta_I)[0];  // we're sure that there is only 1 element in the std::vector<size_t>

                // Calculate the total sign
                int sign = beta_I.operatorPhaseFactor(p) * beta_J.operatorPhaseFactor(q);

                double value = hamiltonian_parameters.get_h()(p,q);

                method(I, J, sign*value);
                method(J, I, sign*value);

                for (size_t r = 0; r < K; r++) {  // r loops over spatial orbitals

                    if (beta_I.isOccupied(r) && beta_J.isOccupied(r)) {  // r must be occupied on the left and on the right
                        if ((p != r) && (q != r)) {  // can't create or annihilate the same orbital
                            double value = hamiltonian_parameters.get_g()(p,q,r,r)
                                           - hamiltonian_parameters.get_g()(r,q,p,r)
                                           - hamiltonian_parameters.get_g()(p,r,r,q)
                                           + hamiltonian_parameters.get_g()(r,r,p,q);

                            method(I, J, sign*value);
                            method(J, I, sign*value);
                        }
                    }

                    if (alpha_I.isOccupied(r)) {  // alpha_I == alpha_J from the previous if-branch

                        double value = hamiltonian_parameters.get_g()(p,q,r,r)
                                       + hamiltonian_parameters.get_g()(q,p,r,r);

                        method(I, J, sign*value);
                        method(J, I, sign*value);
                    }
                }
            }

            // 1 electron excitation in alpha, 1 in beta
            if ((alpha_I.countNumberOfDifferences(alpha_J) == 2) && (beta_I.countNumberOfDifferences(beta_J) == 2)) {

                // Find the orbitals that are occupied in one string, and aren't in the other
                size_t p = alpha_I.findDifferentOccupations(alpha_J)[0];  // we're sure that there is only 1 element in the std::vector<size_t>
                size_t q = alpha_J.findDifferentOccupations(alpha_I)[0];  // we're sure that there is only 1 element in the std::vector<size_t>

                size_t r = beta_I.findDifferentOccupations(beta_J)[0];  // we're sure that there is only 1 element in the std::vector<size_t>
                size_t s = beta_J.findDifferentOccupations(beta_I)[0];  // we're sure that there is only 1 element in the std::vector<size_t>

                int sign = alpha_I.operatorPhaseFactor(p) * alpha_J.operatorPhaseFactor(q) * beta_I.operatorPhaseFactor(r) * beta_J.operatorPhaseFactor(s);
                double value = hamiltonian_parameters.get_g()(p,q,r,s)
                               + hamiltonian_parameters.get_g()(q,p,s,r);

                method(I, J, sign*value);
                method(J, I, sign*value);
            }

            // 2 electron excitations in alpha, 0 in beta
            if ((alpha_I.countNumberOfDifferences(alpha_J) == 4) && (beta_I.countNumberOfDifferences(beta_J) == 0)) {

                // Find the orbitals that are occupied in one string, and aren't in the other
                std::vector<size_t> occupied_indices_I = alpha_I.findDifferentOccupations(alpha_J);  // we're sure this has two elements
                size_t p = occupied_indices_I[0];
                size_t r = occupied_indices_I[1];

                std::vector<size_t> occupied_indices_J = alpha_J.findDifferentOccupations(alpha_I);  // we're sure this has two elements
                size_t q = occupied_indices_J[0];
                size_t s = occupied_indices_J[1];

                int sign = alpha_I.operatorPhaseFactor(p) * alpha_I.operatorPhaseFactor(r) * alpha_J.operatorPhaseFactor(q) * alpha_J.operatorPhaseFactor(s);

                double value = hamiltonian_parameters.get_g()(p,q,r,s)
                            - hamiltonian_parameters.get_g()(p,s,r,q)
                            - hamiltonian_parameters.get_g()(r,q,p,s)
                            + hamiltonian_parameters.get_g()(r,s,p,q);

                method(I, J, sign*value);
                method(J, I, sign*value);
            }

            // 0 electron excitations in alpha, 2 in beta
            if ((alpha_I.countNumberOfDifferences(alpha_J) == 0) && (beta_I.countNumberOfDifferences(beta_J) == 4)) {

                // Find the orbitals that are occupied in one string, and aren't in the other
                std::vector<size_t> occupied_indices_I = beta_I.findDifferentOccupations(beta_J);  // we're sure this has two elements
                size_t p = occupied_indices_I[0];
                size_t r = occupied_indices_I[1];

                std::vector<size_t> occupied_indices_J = beta_J.findDifferentOccupations(beta_I);  // we're sure this has two elements
                size_t q = occupied_indices_J[0];
                size_t s = occupied_indices_J[1];

                int sign = beta_I.operatorPhaseFactor(p) * beta_I.operatorPhaseFactor(r) * beta_J.operatorPhaseFactor(q) * beta_J.operatorPhaseFactor(s);

                double value = hamiltonian_parameters.get_g()(p,q,r,s)
                               - hamiltonian_parameters.get_g()(p,s,r,q)
                               - hamiltonian_parameters.get_g()(r,q,p,s)
                               + hamiltonian_parameters.get_g()(r,s,p,q);

                method(I, J, sign*value);
                method(J, I, sign*value);
            }
        }  // loop over addresses J > I
    }  // loop over addresses I
}


/*
 *  CONSTRUCTORS
 */

/**
 *  @param fock_space       the full alpha and beta product Fock space
 */
SelectedCI::SelectedCI(const ProductFockSpace& fock_space) :
    HamiltonianBuilder(),
    fock_space(fock_space)
{}



/*
 *  OVERRIDDEN PUBLIC METHODS
 */

/**
 *  @param hamiltonian_parameters       the SelectedCI Hamiltonian parameters in an orthonormal orbital basis
 *
 *  @return the SelectedCI Hamiltonian matrix
 */
Eigen::MatrixXd SelectedCI::constructHamiltonian(const HamiltonianParameters& hamiltonian_parameters) const {
    auto K = hamiltonian_parameters.get_h().get_dim();
    if (K != this->fock_space.get_K()) {
        throw std::invalid_argument("Basis functions of the Fock space and hamiltonian_parameters are incompatible.");
    }

    auto dim = fock_space.get_dimension();

    Eigen::MatrixXd result_matrix = Eigen::MatrixXd::Zero(dim, dim);
    result_matrix += this->calculateDiagonal(hamiltonian_parameters).asDiagonal();

    // We pass to a matrix and create the corresponding lambda function
    PassToMethod addToMatrix = [&result_matrix](size_t I, size_t J, double value) { result_matrix(I, J) += value; };

    this->evaluateHamiltonian(hamiltonian_parameters, addToMatrix);
    return result_matrix;
}


/**
 *  @param hamiltonian_parameters       the SelectedCI Hamiltonian parameters in an orthonormal orbital basis
 *  @param x                            the vector upon which the SelectedCI Hamiltonian acts
 *  @param diagonal                     the diagonal of the SelectedCI Hamiltonian matrix
 *
 *  @return the action of the SelectedCI Hamiltonian on the coefficient vector
 */
Eigen::VectorXd SelectedCI::matrixVectorProduct(const HamiltonianParameters& hamiltonian_parameters, const Eigen::VectorXd& x, const Eigen::VectorXd& diagonal) const {

    auto K = hamiltonian_parameters.get_h().get_dim();
    if (K != this->fock_space.get_K()) {
        throw std::invalid_argument("Basis functions of the Fock space and hamiltonian_parameters are incompatible.");
    }

    Eigen::VectorXd matvec = diagonal.cwiseProduct(x);

    // We pass to a the matvec and create the corresponding lambda function
    PassToMethod addToMatvec = [&matvec, &x](size_t I, size_t J, double value) { matvec(I) += value * x(J); };

    this->evaluateHamiltonian(hamiltonian_parameters, addToMatvec);

    return matvec;
}


/**
 *  @param hamiltonian_parameters       the SelectedCI Hamiltonian parameters in an orthonormal orbital basis
 *
 *  @return the diagonal of the matrix representation of the SelectedCI Hamiltonian
 */
Eigen::VectorXd SelectedCI::calculateDiagonal(const HamiltonianParameters &hamiltonian_parameters) const {

    auto K = hamiltonian_parameters.get_h().get_dim();
    if (K != this->fock_space.get_K()) {
        throw std::invalid_argument("Basis functions of the Fock space and hamiltonian_parameters are incompatible.");
    }

    FockSpace fock_space_alpha = fock_space.get_fock_space_alpha();
    FockSpace fock_space_beta = fock_space.get_fock_space_beta();

    auto dim = fock_space.get_dimension();

    // Diagonal contributions
    Eigen::VectorXd diagonal = Eigen::VectorXd::Zero(dim);

    for (size_t I = 0; I < dim; I++) {  // Ia loops over addresses of alpha onvs
        Configuration configuration_I = this->fock_space.get_configuration(I);
        ONV alpha_I = configuration_I.onv_alpha;
        ONV beta_I = configuration_I.onv_beta;

        for (size_t p = 0; p < K; p++) {

            if (alpha_I.isOccupied(p)) {
                diagonal(I) += hamiltonian_parameters.get_h()(p,p);
                for (size_t q = 0; q < K; q++) {
                    if (beta_I.isOccupied(q)) {

                        diagonal(I) += hamiltonian_parameters.get_g()(p,p,q,q);
                    }

                    if (p != q) {  // can't create/annihilate the same orbital twice
                        if (alpha_I.isOccupied(q)) {
                            diagonal(I) += hamiltonian_parameters.get_g()(p,p,q,q);
                            diagonal(I) += hamiltonian_parameters.get_g()(p,p,q,q);
                        }
                    }

                }  // loop over q
            }

            if (beta_I.isOccupied(p)) {
                diagonal(I) += hamiltonian_parameters.get_h()(p,p);
                for (size_t q = 0; q < K; q++) {
                    if (alpha_I.isOccupied(q)) {
                        diagonal(I) += hamiltonian_parameters.get_g()(p,p,q,q);
                    }

                    if (p != q) {  // can't create/annihilate the same orbital twice
                        if (beta_I.isOccupied(q)) {
                            diagonal(I) += hamiltonian_parameters.get_g()(p,p,q,q);
                            diagonal(I) += hamiltonian_parameters.get_g()(p,p,q,q);
                        }
                    }
                }  // loop over q
            }
        }  // loop over q

    }  // alpha address (Ia) loop
    return diagonal;

}



}  // namespace GQCP
