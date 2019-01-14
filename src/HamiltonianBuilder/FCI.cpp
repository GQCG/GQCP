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
#include "HamiltonianBuilder/FCI.hpp"
#include <chrono>

namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param fock_space       the full alpha and beta product Fock space
 */
FCI::FCI(const ProductFockSpace& fock_space) :
        HamiltonianBuilder(),
        fock_space (fock_space)
{
    this->alpha_resolved = this->alphaOneElectronCouplings();
}

/*
 *
 */

/**
 *  Calculates all Hamiltonian elements for operators exclusively operating for one spin function
 *  and stores these in a sparse matrix
 *
 *  @param fock_space                   Fock space for the spin function specific Hamiltonian
 *  @param k                            Modified one-electron operator
 *  @param hamiltonian_parameters       The Hamiltonian parameters in an orthonormal orbital basis
 *  @param sparse_mat                   The representation of the spin function specific Hamiltonian
 */
void FCI::spinSeparatedModule(FockSpace& fock_space, const OneElectronOperator& k,
                              const HamiltonianParameters& hamiltonian_parameters, Eigen::SparseMatrix<double>& sparse_mat){

    size_t K = fock_space.get_K();
    size_t N = fock_space.get_N();
    size_t dim = fock_space.get_dimension();
    
    std::vector<Eigen::Triplet<double>> triplet_vector;
    triplet_vector.reserve(fock_space.totalTwoElectronCouplingCount());
    
    ONV onv = fock_space.get_ONV(0);  // onv with address 0
    for (size_t I = 0; I < dim; I++) {  // I loops over all addresses in the Fock space
        if (I > 0) {
            fock_space.setNext(onv);
        }
        int sign1 = -1;
        for (size_t e1 = 0; e1 < N; e1++) {  // A1 (annihilation 1)

            sign1 *= -1;
            size_t p = onv.get_occupied_index(e1);  // retrieve the index of a given electron
            size_t address = I - fock_space.get_vertex_weights(p, e1 + 1);


            size_t address1 = address;
            size_t e2 = e1;
            size_t q = p;
            /**
             *  A1 > C1 (annihlation 1 > creation 1)
             */
            int sign2 = sign1;
            q--;
            e2--;
            fock_space.shiftUntilPreviousUnoccupiedOrbital<1>(onv, address1, q, e2, sign2);
            while (q != -1) {

                size_t address2 = address1 + fock_space.get_vertex_weights(q, e2 + 2);

                /**
                 *  A2 > C2
                 */
                int sign3 = sign1;
                for (size_t e3 = e1 + 1; e3 < N; e3++) {
                    sign3 *= -1;  // initial sign3 = sign of the annhilation, with one extra electron(from crea) = *-1
                    size_t r = onv.get_occupied_index(e3);
                    size_t address3 = address2 - fock_space.get_vertex_weights(r, e3 + 1);

                    size_t e4 = e3 + 1;
                    size_t s = r + 1;

                    int sign4 = sign3;
                    fock_space.shiftUntilNextUnoccupiedOrbital<1>(onv, address3, s, e4, sign4);

                    while (s < K) {
                        size_t J = address3 + fock_space.get_vertex_weights(s, e4);
                        int signev = sign1 * sign2 * sign3 * sign4;
                        double value = signev * 0.5 * (hamiltonian_parameters.get_g()(p, q, r, s) +
                                                       hamiltonian_parameters.get_g()(r, s, p, q) -
                                                       hamiltonian_parameters.get_g()(p, s, r, q) -
                                                       hamiltonian_parameters.get_g()(r, q, p, s));


                        triplet_vector.emplace_back(I,J, value);
                        triplet_vector.emplace_back(J,I, value);

                        s++;
                        fock_space.shiftUntilNextUnoccupiedOrbital<1>(onv, address3, s, e4, sign4);
                    }
                }

                q--;
                fock_space.shiftUntilPreviousUnoccupiedOrbital<1>(onv, address1, q, e2, sign2);

            }



            e2 = e1 + 1;
            q = p + 1;
            sign2 = sign1;
            fock_space.shiftUntilNextUnoccupiedOrbital<1>(onv, address, q, e2, sign2);


            /**
             *  A1 < C1
             */
            while (q < K) {

                // BRANCH N


                address1 = address + fock_space.get_vertex_weights(q, e2);
                /**
                 *  A2 > C1
                 */
                int sign3 = sign2;
                for (size_t e3 = e2; e3 < N; e3++) {
                    sign3 *= -1; // -1 cause we created electron (creation) sign of A is now the that of C *-1
                    size_t r = onv.get_occupied_index(e3);
                    size_t address3 = address1 - fock_space.get_vertex_weights(r, e3 + 1);

                    size_t e4 = e3 + 1;
                    size_t s = r + 1;

                    
                    int sign4 = sign3;
                    fock_space.shiftUntilNextUnoccupiedOrbital<1>(onv, address3, s, e4, sign4);

                    while (s < K) {
                        size_t J = address3 + fock_space.get_vertex_weights(s, e4);
                        int signev = sign1 * sign2 * sign3 * sign4;

                        double value = signev * 0.5 * (hamiltonian_parameters.get_g()(p, q, r, s) +
                                                       hamiltonian_parameters.get_g()(r, s, p, q) -
                                                       hamiltonian_parameters.get_g()(r, q, p, s) -
                                                       hamiltonian_parameters.get_g()(p, s, r, q));

                        triplet_vector.emplace_back(I,J, value);
                        triplet_vector.emplace_back(J,I, value);

                        s++;  // go to the next orbital
                        fock_space.shiftUntilNextUnoccupiedOrbital<1>(onv, address3, s, e4, sign4);

                    }  // (creation)

                }


                size_t r = q;

                sign3 = sign2;

                size_t address1c = address1;

                /**
                 *  A2 < C1, (A2 > A1)
                 */
                for (size_t e3 = e2 - 1; e3 > e1; e3--) {
                    sign3 *= -1;
                    size_t e4 = e2;
                    address1c += fock_space.get_vertex_weights(r, e3) -
                                 fock_space.get_vertex_weights(r, e3 + 1);
                    r = onv.get_occupied_index(e3);
                    size_t address2 = address1c - fock_space.get_vertex_weights(r, e3);
                    int sign4 = sign2;
                    size_t s = q + 1;
                    fock_space.shiftUntilNextUnoccupiedOrbital<1>(onv, address2, s, e4, sign4);
                    while (s < K) {

                        size_t J = address2 + fock_space.get_vertex_weights(s, e4);

                        int signev = sign1 * sign2 * sign3 * sign4;
                        double value = signev * 0.5 * (hamiltonian_parameters.get_g()(p, q, r, s) +
                                                       hamiltonian_parameters.get_g()(r, s, p, q) -
                                                       hamiltonian_parameters.get_g()(r, q, p, s) -
                                                       hamiltonian_parameters.get_g()(p, s, r, q));

                        triplet_vector.emplace_back(I,J, value);
                        triplet_vector.emplace_back(J,I, value);

                        s++;

                        fock_space.shiftUntilNextUnoccupiedOrbital<1>(onv, address2, s, e4, sign4);

                    }


                }

                /**
                 *  A2 = C1
                 */
                int signev = sign2 * sign1;
                
                double value_I =  k(p, q);  // cover the one electron calculations
                
                for (size_t s = 0; s < K; s++) {
                    if(!onv.isOccupied(s)){

                        value_I += 0.5 * (hamiltonian_parameters.get_g()(p, s, s, q));

                    } else {

                        value_I  += 0.5 *  (hamiltonian_parameters.get_g()(s, s, p, q) - hamiltonian_parameters.get_g()(s, q, p, s) + hamiltonian_parameters.get_g()(p, q, s, s));

                    }

                }

                value_I *= signev;


                q++;



                triplet_vector.emplace_back(I,address1, value_I);
                triplet_vector.emplace_back(address1,I, value_I);

                fock_space.shiftUntilNextUnoccupiedOrbital<1>(onv, address, q, e2, sign2);


            }
        }
    }

    sparse_mat.setFromTriplets(triplet_vector.begin(),triplet_vector.end());


}

/**
 *  Calculates all one-electron couplings for the beta Fock space
 *  and attributes two-electron integrals based on the one-electron indexes of the coupling and two fixed indexes
 *
 *  @param r                        Fixed index of two-electorn integral
 *  @param s                        Fixed index of two-electron integral
 *  @param hamiltonian_parameters   The Hamiltonian parameters in an orthonormal orbital basis
 *
 *  @return                         The sparse matrix containing the calculated two-electron integrals mapped to one-electron couplings
 */
Eigen::SparseMatrix<double> FCI::betaTwoElectronOneElectronModule(size_t r, size_t s,
                                                                  const HamiltonianParameters& hamiltonian_parameters) {

    FockSpace beta = fock_space.get_fock_space_beta();
    const bool do_diagonal = (r != s);
    size_t K = beta.get_K();
    size_t N = beta.get_N();
    size_t dim = beta.get_dimension();
    Eigen::SparseMatrix<double> sparseMatrix(dim, dim);
    std::vector<Eigen::Triplet<double>> triplet_vector;
    size_t mod = 0;
    if (do_diagonal){
        mod += dim;
    }
    triplet_vector.reserve(beta.totalOneElectronCouplingCount() + mod);
    ONV onv = beta.get_ONV(0);  // onv with address 0
    for (size_t I = 0; I < dim; I++) {  // I loops over all the addresses of the onv
        for (size_t e1 = 0; e1 < N; e1++) {  // e1 (electron 1) loops over the (number of) electrons
            size_t p = onv.get_occupied_index(e1);  // retrieve the index of a given electron
            // remove the weight from the initial address I, because we annihilate
            size_t address = I - beta.get_vertex_weights(p, e1 + 1);

            if(do_diagonal){
                triplet_vector.emplace_back(I, I, hamiltonian_parameters.get_g()(r, s, p, p));
            }

            // The e2 iteration counts the amount of encountered electrons for the creation operator
            // We only consider greater addresses than the initial one (because of symmetry)
            // Hence we only count electron after the annihilated electron (e1)
            size_t e2 = e1 + 1;
            size_t q = p + 1;

            int sign_e2 = 1;
            // perform a shift
            beta.shiftUntilNextUnoccupiedOrbital<1>(onv, address, q, e2, sign_e2);

            while (q < K) {
                size_t J = address + beta.get_vertex_weights(q, e2);
                double value = sign_e2*hamiltonian_parameters.get_g()(r, s, p, q);
                triplet_vector.emplace_back(I, J, value);
                triplet_vector.emplace_back(J, I, value);

                q++; // go to the next orbital

                // perform a shift
                beta.shiftUntilNextUnoccupiedOrbital<1>(onv, address, q, e2, sign_e2);
            }  //  (creation)


        } // e1 loop (annihilation)


        // Prevent last permutation
        if (I < dim - 1) {
            beta.setNext(onv);
        }
    }
    sparseMatrix.setFromTriplets(triplet_vector.begin(),triplet_vector.end());
    return sparseMatrix;
    
}



/**
 *  Calculates all one-eletron couplings for each annihilation-creation pair in the alpha Fock space
 *  and stores them in sparse matrices for each combination
 *
 *  @return vector of sparse matrices containing the one-electron couplings for the alpha Fock space
 */
std::vector<Eigen::SparseMatrix<double>> FCI::alphaOneElectronCouplings() {

    GQCP::FockSpace alpha = this->fock_space.get_fock_space_alpha();

    size_t K = alpha.get_K();
    size_t N = alpha.get_N();
    size_t dim = alpha.get_dimension();

    std::vector<std::vector<Eigen::Triplet<double>>> sparseEntries(K*(K+1)/2);
    std::vector<Eigen::SparseMatrix<double>> matrixes(K*(K+1)/2, Eigen::SparseMatrix<double>(dim, dim));
    for (std::vector<Eigen::Triplet<double>>& reserver : sparseEntries) {
        reserver.reserve(2*(1+N*(K-N)));
    }


    ONV onv = alpha.get_ONV(0);  // onv with address 0
    for (size_t I = 0; I < dim; I++) {  // I loops over all the addresses of the onv
        for (size_t e1 = 0; e1 < N; e1++) {  // e1 (electron 1) loops over the (number of) electrons
            size_t p = onv.get_occupied_index(e1);  // retrieve the index of a given electron
            // remove the weight from the initial address I, because we annihilate
            size_t address = I - alpha.get_vertex_weights(p, e1 + 1);
            // The e2 iteration counts the amount of encountered electrons for the creation operator
            // We only consider greater addresses than the initial one (because of symmetry)
            // Hence we only count electron after the annihilated electron (e1)
            sparseEntries[p*(K+K+1-p)/2].emplace_back(I, I, 1);
            size_t e2 = e1 + 1;
            size_t q = p + 1;
            int sign_e2 = 1;
            // perform a shift
            alpha.shiftUntilNextUnoccupiedOrbital<1>(onv, address, q, e2, sign_e2);
            size_t dex = 0;
            while (q < K) {
                size_t J = address + alpha.get_vertex_weights(q, e2);

                sparseEntries[p*(K+K+1-p)/2 + q - p].emplace_back(I, J, sign_e2);
                sparseEntries[p*(K+K+1-p)/2 + q - p].emplace_back(J, I, sign_e2);

                q++; // go to the next orbital
                dex++;
                // perform a shift
                alpha.shiftUntilNextUnoccupiedOrbital<1>(onv, address, q, e2, sign_e2);
            }  //  (creation)


        } // e1 loop (annihilation)


        // Prevent last permutation
        if (I < dim - 1) {
            alpha.setNext(onv);
        }
    }


    for (size_t k = 0; k < K*(K+1)/2 ; k++){
        matrixes[k].setFromTriplets(sparseEntries[k].begin(), sparseEntries[k].end());
    }

    return matrixes;


}



/*
 *  OVERRIDDEN PUBLIC METHODS
 */

/**
 *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
 *
 *  @return the FCI Hamiltonian matrix
 */
Eigen::MatrixXd FCI::constructHamiltonian(const HamiltonianParameters& hamiltonian_parameters) const {

    auto K = hamiltonian_parameters.get_h().get_dim();
    if (K != this->fock_space.get_K()) {
        throw std::invalid_argument("Basis functions of the Fock space and hamiltonian_parameters are incompatible.");
    }

    Eigen::MatrixXd result_matrix = Eigen::MatrixXd::Zero(this->fock_space.get_dimension(), this->fock_space.get_dimension());

    FockSpace fock_space_alpha = fock_space.get_fock_space_alpha();
    FockSpace fock_space_beta = fock_space.get_fock_space_beta();

    OneElectronOperator k = hamiltonian_parameters.calculateEffectiveOneElectronIntegrals();

    auto N_alpha = fock_space_alpha.get_N();
    auto dim_alpha = fock_space_alpha.get_dimension();
    auto N_beta = fock_space_beta.get_N();
    auto dim_beta = fock_space_beta.get_dimension();

    if (!is_ham_par_set) {
        spinSeparatedModule(fock_space_alpha, k, hamiltonian_parameters, alpha_ev);
        spinSeparatedModule(fock_space_beta, k, hamiltonian_parameters, beta_ev);
        for (size_t p = 0; p < K; p++) {
            beta_resolved[p * (K + K + 1 - p) / 2] = betaTwoElectronOneElectronModule(p, p, hamiltonian_parameters);
            for (size_t q = p + 1; q < K; q++) {
                beta_resolved[p * (K + K + 1 - p) / 2 + q - p] = betaTwoElectronOneElectronModule(p, q, hamiltonian_parameters);
            }
        }

    }

    for (size_t i = 0; i < dim_alpha; i++) {
        for (size_t j = 0; j < dim_alpha; j++) {
            result_matrix.block(i * dim_beta, j * dim_beta, dim_beta, dim_beta) += beta_ev;
        }
    }

    Eigen::MatrixXd ones = Eigen::MatrixXd::Ones(dim_alpha, dim_alpha);

    for (int i = 0; i < alpha_ev.outerSize(); ++i){
        for (Eigen::SparseMatrix<double>::InnerIterator it(alpha_ev, i); it; ++it) {
            result_matrix.block(it.row() * dim_alpha, it.col() * dim_alpha, dim_alpha, dim_alhpa) += it_value*ones;
        }
    }


    for (const Eigen::SparseMatrix<double> alpha_couplings : alpha_resolved) {
        for (i = 0; i < alpha_couplings.outerSize(); ++i){
            for (Eigen::SparseMatrix<double>::InnerIterator it(alpha_couplings, i); it; ++it) {
                result_matrix.block(it.row() * dim_alpha, it.col() * dim_alpha, dim_alpha, dim_alpha) += it_value*ones;
            }
        }
    }


    for (size_t p = 0; p<K; p++) {
        matvecmap.noalias() += alpha_resolved[p*(K+K+1-p)/2] * xmap * beta_resolved[p*(K+K+1-p)/2];
        for (size_t q = p + 1; q<K; q++) {
            matvecmap.noalias() += alpha_resolved[p*(K+K+1-p)/2 + q - p] * xmap * beta_resolved[p*(K+K+1-p)/2 + q - p];
        }
    }

    return result_matrix;
}


/**
 *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
 *  @param x                            the vector upon which the FCI Hamiltonian acts
 *  @param diagonal                     the diagonal of the FCI Hamiltonian matrix
 *
 *  @return the action of the FCI Hamiltonian on the coefficient vector
 */
Eigen::VectorXd FCI::matrixVectorProduct(const HamiltonianParameters& hamiltonian_parameters, const Eigen::VectorXd& x, const Eigen::VectorXd& diagonal) const {
    auto K = hamiltonian_parameters.get_h().get_dim();
    if (K != this->fock_space.get_K()) {
        throw std::invalid_argument("Basis functions of the Fock space and hamiltonian_parameters are incompatible.");
    }

    FockSpace fock_space_alpha = fock_space.get_fock_space_alpha();
    FockSpace fock_space_beta = fock_space.get_fock_space_beta();

    auto dim_alpha = fock_space_alpha.get_dimension();
    auto dim_beta = fock_space_beta.get_dimension();

    Eigen::VectorXd matvec = diagonal.cwiseProduct(x);

    // Calculate the effective one-electron integrals
    OneElectronOperator k = hamiltonian_parameters.calculateEffectiveOneElectronIntegrals();

    GQCP::OneElectronOperator k = hamiltonian_parameters.calculateEffectiveOneElectronIntegrals();
    Eigen::Map<Eigen::MatrixXd> matvecmap(matvec.data(), dim_alpha, dim_beta);
    Eigen::Map<const Eigen::MatrixXd> xmap(x.data(), dim_alpha, dim_beta);




    if (xyr) {
        xyr = false;

        this->spinSeparatedModule(fock_space_alpha, k, hamiltonian_parameters, alpha_ev);
        this->spinSeparatedModule(fock_space_beta, k, hamiltonian_parameters, beta_ev);

        beta_resolved = std::vector<Eigen::SparseMatrix<double>>(K * (K + 1) / 2,
                                                                 Eigen::SparseMatrix<double>(dim_beta, dim_beta));

        for (size_t p = 0; p < K; p++) {
            beta_resolved[p * (K + K + 1 - p) / 2] = betaTwoElectronOneElectronModule(p, p, hamiltonian_parameters);
            for (size_t q = p + 1; q < K; q++) {
                beta_resolved[p * (K + K + 1 - p) / 2 + q - p] = betaTwoElectronOneElectronModule(p, q,
                                                                                                  hamiltonian_parameters);
            }
        }
    }



    for (size_t p = 0; p<K; p++) {
        matvecmap.noalias() += alpha_resolved[p*(K+K+1-p)/2] * xmap * beta_resolved[p*(K+K+1-p)/2];
        for (size_t q = p + 1; q<K; q++) {
            matvecmap.noalias() += alpha_resolved[p*(K+K+1-p)/2 + q - p] * xmap * beta_resolved[p*(K+K+1-p)/2 + q - p];
        }
    }



    matvecmap.noalias() += alpha_ev * xmap + xmap * beta_ev;

    return matvec;
}


/**
 *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
 *
 *  @return the diagonal of the matrix representation of the Hamiltonian
 */
Eigen::VectorXd FCI::calculateDiagonal(const HamiltonianParameters& hamiltonian_parameters) const {

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
    Eigen::VectorXd diagonal =  Eigen::VectorXd::Zero(dim);

    // Calculate the effective one-electron integrals
    OneElectronOperator k = hamiltonian_parameters.calculateEffectiveOneElectronIntegrals();

    ONV spin_string_alpha = fock_space_alpha.makeONV(0);
    for (size_t Ia = 0; Ia < dim_alpha; Ia++) {  // Ia loops over addresses of alpha spin strings

        ONV spin_string_beta = fock_space_beta.makeONV(0);
        for (size_t Ib = 0; Ib < dim_beta; Ib++) {  // Ib loops over addresses of beta spin strings

            for (size_t p = 0; p < K; p++) {  // p loops over SOs

                if (spin_string_alpha.isOccupied(p)) {  // p is in Ia
                    diagonal(Ia * dim_beta + Ib) += k(p, p);

                    for (size_t q = 0; q < K; q++) {  // q loops over SOs
                        if (spin_string_alpha.isOccupied(q)) {  // q is in Ia
                            diagonal(Ia * dim_beta + Ib) += 0.5 * hamiltonian_parameters.get_g()(p, p, q, q);
                        } else {  // q is not in I_alpha
                            diagonal(Ia * dim_beta + Ib) += 0.5 * hamiltonian_parameters.get_g()(p, q, q, p);
                        }

                        if (spin_string_beta.isOccupied(q)) {  // q is in Ib
                            diagonal(Ia * dim_beta + Ib) += hamiltonian_parameters.get_g()(p, p, q, q);
                        }
                    }  // q loop
                }


                if (spin_string_beta.isOccupied(p)) {  // p is in Ib
                    diagonal(Ia * dim_beta + Ib) += k(p, p);


                    for (size_t q = 0; q < K; q++) {  // q loops over SOs
                        if (spin_string_beta.isOccupied(q)) {  // q is in Ib
                            diagonal(Ia * dim_beta + Ib) += 0.5 * hamiltonian_parameters.get_g()(p, p, q, q);

                        } else {  // q is not in I_beta
                            diagonal(Ia * dim_beta + Ib) += 0.5 * hamiltonian_parameters.get_g()(p, q, q, p);
                        }
                    }  // q loop
                }

            }  // p loop

            if (Ib < dim_beta - 1) {  // prevent last permutation to occur
                fock_space_beta.setNextONV(spin_string_beta);
            }
        }  // beta address (Ib) loop

        if (Ia < dim_alpha - 1) {  // prevent last permutation to occur
            fock_space_alpha.setNextONV(spin_string_alpha);
        }
    }  // alpha address (Ia) loop

    return diagonal;
}


}  // namespace GQCP
