#include "HamiltonianBuilder/FCI.hpp"


namespace GQCG {


/*
 *  CONSTRUCTORS
 */

/**
 *  Constructor given a @param hamiltonian_parameters and @param fock_space
 */
FCI::FCI(FockSpaceProduct fock_space) :
        HamiltonianBuilder (),
        fock_space (fock_space)
{}



/*
 *  OVERRIDDEN PUBLIC METHODS
 */

/**
 *  @return Hamiltonian matrix as an Eigen::MatrixXd given @param hamiltonian_parameters
 */
Eigen::MatrixXd FCI::constructHamiltonian(const HamiltonianParameters& hamiltonian_parameters) {

    auto K = hamiltonian_parameters.get_h().get_dim();
    if (K != this->fock_space.get_K()) {
        throw std::invalid_argument("Basis functions of the Fock space and hamiltonian_parameters are incompatible.");
    }

    Eigen::MatrixXd result_matrix = Eigen::MatrixXd::Zero(this->fock_space.get_dimension(), this->fock_space.get_dimension());
    
    FockSpace fock_space_alpha = fock_space.get_fock_space_alpha();
    FockSpace fock_space_beta = fock_space.get_fock_space_beta();

    auto N_alpha = fock_space_alpha.get_N();
    auto dim_alpha = fock_space_alpha.get_dimension();
    auto N_beta = fock_space_beta.get_N();
    auto dim_beta = fock_space_beta.get_dimension();

    alpha_one_electron_couplings = { dim_alpha, std::vector<OneElectronCoupling>(N_alpha * (K + 1 - N_alpha)) };
    beta_one_electron_couplings = { dim_beta, std::vector<OneElectronCoupling>(N_beta * (K + 1 - N_beta)) };

    // 1. ALPHA-ALPHA
    ONV spin_string_alpha = fock_space_alpha.get_ONV(0);  // alpha spin string with address 0
    for (size_t I_alpha = 0; I_alpha < dim_alpha; I_alpha++) {  // I_alpha loops over all the addresses of the alpha spin strings

        size_t coupling_address_index = 0;  // index of |J_alpha> in the (N_alpha * (K + 1 - N_alpha))-long std::vector
        // located at alpha_one_electron_couplings[I_alpha]

        for (size_t p = 0; p < K; p++) {  // p loops over SOs
            int sign_p = 1;  // sign for the annihilation operator (a_p)

            if (spin_string_alpha.annihilate(p, sign_p)) {
                for (size_t q = 0; q < K; q++) {  // q loops over SOs

                    // one-electron contributions for alpha, i.e. one electron excitation
                    int sign_pq = sign_p;  // sign for the total excitation operator (a^\dagger_q a_p)
                    if (spin_string_alpha.create(q, sign_pq)) {

                        size_t J_alpha = fock_space_alpha.getAddress(spin_string_alpha);

                        // For the 'diagonal beta contributions', i.e. I_beta = J_beta, the one-electron alpha contributions
                        // are the same
                        // We are storing the alpha addresses as 'major', i.e. the total address I_alpha I_beta = I_alpha * dim_beta + I_beta
                        for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {
                            double value = sign_pq * hamiltonian_parameters.get_h().get(p,q);
                            result_matrix(I_alpha * dim_beta + I_beta, J_alpha * dim_beta + I_beta) += value;
                        }

                        // We have found a spin string that is one electron excitation away from |I_alpha>
                        // We will store it, since these strings are also needed in the alpha-beta part
                        this->alpha_one_electron_couplings[I_alpha][coupling_address_index] = OneElectronCoupling{sign_pq, p, q, J_alpha};
                        coupling_address_index++;
                        spin_string_alpha.annihilate(q);  // undo the previous creation on q
                    }  // create on q (alpha)


                    // two-electron contributions for beta-beta, i.e. two electron excitations
                    sign_pq = sign_p;  // sign for the total excitation operator (a^\dagger_q a_p)
                    // we have to reset this because we changed this in the previous if-statement

                    if (spin_string_alpha.annihilate(q, sign_pq)) {

                        for (size_t r = 0; r < K; r++) {
                            int sign_pqr = sign_pq;  // sign for total operator (a^\dagger_r a_q a_p)

                            if (spin_string_alpha.create(r, sign_pqr)) {
                                for (size_t s = 0; s < K; s++) {

                                    int sign_pqrs = sign_pqr;  // sign for total operator (a^dagger_s a^\dagger_r a_q a_p)
                                    if (spin_string_alpha.create(s, sign_pqrs)) {

                                        size_t Ja = fock_space_alpha.getAddress(spin_string_alpha);

                                        // For the 'diagonal beta contributions', i.e. Ib = Jb, the two-electron alpha
                                        // contributions are the same

                                        // We are storing the alpha addresses as 'major', i.e. the total address I_alpha I_beta = I_alpha * dim_b + I_b
                                        for (size_t Ib = 0; Ib < dim_beta; Ib++) {
                                            double value = sign_pqrs * 0.5 * hamiltonian_parameters.get_g().get(s, p, r, q);
                                            result_matrix(I_alpha * dim_beta + Ib, Ja * dim_beta + Ib) += value;
                                        }

                                        spin_string_alpha.annihilate(s);  // undo the previous creation on s
                                    }  // create on s (alpha)
                                }  // loop over s

                                spin_string_alpha.annihilate(r);  // undo the previous creation on r
                            }  // create on r (alpha)
                        }  // loop over r

                        spin_string_alpha.create(q);  // undo the previous annihilation on q
                    }  // annihilate on q (alpha)
                }  // loop over q

                spin_string_alpha.create(p);  // undo the previous annihilation on p
            }  // annihilate p (alpha)
        }  // loop over p


        if (I_alpha < dim_alpha - 1) {  // prevent the last permutation to occur
            fock_space_alpha.setNext(spin_string_alpha);
        }
    }  // loop over alpha addresses (I_alpha)


    // 2. BETA-BETA
    ONV spin_string_beta = fock_space_beta.get_ONV(0);  // beta spin string with address 0

    for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {  // I_beta loops over addresses of all beta spin strings

        size_t coupling_address_index = 0;  // index of |J_beta> in the (N_beta * (K + 1 - N_beta))-long std::vector
        // located at alpha_one_electron_couplings[I_alpha]

        for (size_t p = 0; p < K; p++) {  // p loops over SOs
            int sign_p = 1;

            if (spin_string_beta.annihilate(p, sign_p)) {
                for (size_t q = 0; q < K; q++) {  // q loops over SOs

                    // one-electron contributions for beta, i.e. one electron excitation
                    int sign_pq = sign_p;  // sign for the total excitation operator (a^\dagger_q a_p)
                    if (spin_string_beta.create(q, sign_pq)) {

                        size_t J_beta = fock_space_beta.getAddress(spin_string_beta);

                        // For the 'diagonal alpha contributions', i.e. I_alpha = J_alpha, the one-electron beta contributions are
                        // the same
                        // We are storing the alpha addresses as 'major', i.e. the total address I_alpha I_beta = I_alpha * dim_beta + I_beta

                        for (size_t I_alpha = 0; I_alpha < dim_alpha; I_alpha++) {
                            double value = sign_pq * hamiltonian_parameters.get_h().get(p,q);
                            result_matrix (I_alpha * dim_beta + I_beta, I_alpha * dim_beta + J_beta) += value;
                        }

                        // We have found a spin string that is one electron excitation away from |I_alpha>
                        // We will store it, since these strings are also needed in the alpha-beta part
                        this->beta_one_electron_couplings[I_beta][coupling_address_index] = OneElectronCoupling{sign_pq, p, q, J_beta};
                        coupling_address_index++;
                        spin_string_beta.annihilate(q);  // undo the previous creation on q
                    }  // create on q (beta)


                    // two-electron contributions for beta-beta, i.e. two electron excitations
                    sign_pq = sign_p;  // sign for the total excitation operator (a^\dagger_q a_p)
                    // we have to reset this because we changed this in the previous if-statement

                    if (spin_string_beta.annihilate(q, sign_pq)) {

                        for (size_t r = 0; r < K; r++) {  // r loops over SOs
                            int sign_pqr = sign_pq;  // sign for total operator (a^\dagger_r a_q a_p)

                            if (spin_string_beta.create(r, sign_pqr)) {
                                for (size_t s = 0; s < K; s++) {  // s loops over SOs

                                    int sign_pqrs = sign_pqr;  // sign for total operator (a^dagger_s a^\dagger_r a_q a_p)
                                    if (spin_string_beta.create(s, sign_pqrs)) {

                                        size_t Jb = fock_space_beta.getAddress(spin_string_beta);

                                        // For the 'diagonal alpha contributions', i.e. Ia = Ja, the two-electron beta
                                        // contributions are the same

                                        // We are storing the alpha addresses as 'major', i.e. the total address IaIb = Ia * dim_b + I_b
                                        for (size_t Ia = 0; Ia < dim_alpha; Ia++) {
                                            double value = sign_pqrs * 0.5 * hamiltonian_parameters.get_g().get(s, p, r, q);
                                            result_matrix(Ia * dim_beta + I_beta, Ia * dim_beta + Jb) += value;
                                        }

                                        spin_string_beta.annihilate(s);  // undo the previous creation on s
                                    }  // create on s (beta)
                                }  // loop over s

                                spin_string_beta.annihilate(r);  // undo the previous creation on r
                            }  // create on r (beta)
                        }  // loop over r

                        spin_string_beta.create(q);  // undo the previous annihilation on q
                    }  // annihilate on q (beta)
                }  // loop over q

                spin_string_beta.create(p);  // undo the previous annihilation on p
            } // annihilate on p (beta)
        }  // loop over p

        if (I_beta < dim_beta - 1) {  // prevent last permutation to occur
            fock_space_beta.setNext(spin_string_beta);
        }
    }  // loop over beta addresses (I_beta)


    // 3. ALPHA-BETA
    for (size_t I_alpha = 0; I_alpha < dim_alpha; I_alpha++) {  // loop over alpha addresses
        for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {  // loop over beta addresses

            for (const auto& alpha : this->alpha_one_electron_couplings[I_alpha]) {  // traverse all OneElectronCouplings for I_alpha
                for (const auto& beta : this->beta_one_electron_couplings[I_beta]) {  // traverse all OneElectronCouplings for I_beta

                    int sign = alpha.sign * beta.sign;
                    double value = sign * hamiltonian_parameters.get_g().get(alpha.p,alpha.q,beta.p,beta.q);
                    result_matrix( I_alpha * dim_beta + I_beta, alpha.address * dim_beta + beta.address) += value;  // alpha is the major index
                }  // beta OneElectronCouplings
            }  // alpha OneElectronCouplings

        }  // loop over beta addresses (I_beta)
    }  // loop over alpha addresses (I_alpha)
    return result_matrix;
}


/**
 *  @return the action of the Hamiltonian (@param hamiltonian_parameters and @param diagonal) on the coefficient vector @param x
 */
Eigen::VectorXd FCI::matrixVectorProduct(const HamiltonianParameters& hamiltonian_parameters, const Eigen::VectorXd& x, const Eigen::VectorXd& diagonal) {
    auto K = hamiltonian_parameters.get_h().get_dim();
    if (K != this->fock_space.get_K()) {
        throw std::invalid_argument("Basis functions of the Fock space and hamiltonian_parameters are incompatible.");
    }

    FockSpace fock_space_alpha = fock_space.get_fock_space_alpha();
    FockSpace fock_space_beta = fock_space.get_fock_space_beta();

    auto dim_alpha = fock_space_alpha.get_dimension();
    auto dim_beta = fock_space_beta.get_dimension();
    auto dim = fock_space.get_dimension();

    // TODO: use diagonal
    Eigen::VectorXd matvec =  Eigen::VectorXd::Zero(dim);

    // Calculate the effective one-electron integrals
    // TODO: move this to libwint
    Eigen::MatrixXd k_SO = hamiltonian_parameters.get_h().get_matrix_representation();
    for (size_t p = 0; p < K; p++) {
        for (size_t q = 0; q < K; q++) {
            for (size_t r = 0; r < K; r++) {
                k_SO(p,q) -= 0.5 * hamiltonian_parameters.get_g().get(p,r,r,q);
            }
        }
    }


    // ALPHA-ALPHA
    ONV spin_string_alpha_aa = fock_space_alpha.get_ONV(0);  // spin string with address 0
    for (size_t I_alpha = 0; I_alpha < dim_alpha; I_alpha++) {  // I_alpha loops over all the addresses of the alpha spin strings
        if (I_alpha > 0) {
            fock_space_alpha.setNext(spin_string_alpha_aa);
        }

        for (size_t p = 0; p < K; p++) {  // p loops over SOs
            int sign_p = 1;  // sign of the operator a_p_alpha
            if (spin_string_alpha_aa.annihilate(p, sign_p)) {  // if p is in I_alpha

                for (size_t q = 0; q < K; q++) {  // q loops over SOs
                    int sign_pq = sign_p;  // sign of the operator a^dagger_q_alpha a_p_alpha
                    if (spin_string_alpha_aa.create(q, sign_pq)) {  // if q is not occupied in I_alpha
                        size_t J_alpha = fock_space_alpha.getAddress(spin_string_alpha_aa); // find all strings J_alpha that couple to I_alpha

                        for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {  // I_beta loops over all addresses of the beta spin strings
                            matvec(I_alpha*dim_beta + I_beta) += k_SO(p,q) * sign_pq * x(J_alpha*dim_beta + I_beta);  // alpha addresses are major
                        }

                        spin_string_alpha_aa.annihilate(q);  // undo the previous creation
                    }
                }  // q loop

                spin_string_alpha_aa.create(p);  // undo the previous annihilation
            }
        }  // p loop
    }  // I_alpha loop


    // BETA-BETA
    ONV spin_string_beta_bb = fock_space_beta.get_ONV(0);  // spin string with address 0
    for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {  // I_beta loops over all the addresses of the beta spin strings
        if (I_beta > 0) {
            fock_space_beta.setNext(spin_string_beta_bb);
        }

        for (size_t p = 0; p < K; p++) {  // p loops over SOs
            int sign_p = 1;  // sign of the operator a_p_beta
            if (spin_string_beta_bb.annihilate(p, sign_p)) {  // if p is in I_beta

                for (size_t q = 0; q < K; q++) {  // q loops over SOs
                    int sign_pq = sign_p;  // sign of the operator a^dagger_q_beta a_p_beta
                    if (spin_string_beta_bb.create(q, sign_pq)) {  // if q is not occupied in I_beta
                        size_t J_beta = fock_space_beta.getAddress(spin_string_beta_bb);  // find all strings J_beta that couple to I_beta

                        for (size_t I_alpha = 0; I_alpha < dim_alpha; I_alpha++) {  // I_alpha loops over all addresses of the alpha spin strings
                            matvec(I_alpha*dim_beta + I_beta) += k_SO(p,q) * sign_pq * x(I_alpha*dim_beta + J_beta);  // alpha addresses are major
                        }

                        spin_string_beta_bb.annihilate(q);  // undo the previous creation
                    }
                }  // q loop

                spin_string_beta_bb.create(p);  // undo the previous annihilation
            }
        }  // p loop
    }  // I_beta loop

    // ALPHA-ALPHA-ALPHA-ALPHA
    ONV spin_string_alpha_aaaa = fock_space_alpha.get_ONV(0);  // spin string with address 0
    for (size_t I_alpha = 0; I_alpha < dim_alpha; I_alpha++) {  // I_alpha loops over all addresses of alpha spin strings
        if (I_alpha > 0) {
            fock_space_alpha.setNext(spin_string_alpha_aaaa);
        }

        for (size_t p = 0; p < K; p++) {  // p loops over SOs
            int sign_p = 1;  // sign of the operator a_p_alpha
            if (spin_string_alpha_aaaa.annihilate(p, sign_p)) {

                for (size_t q = 0; q < K; q++) {  // q loops over SOs
                    int sign_pq = sign_p;  // sign of the operator a^dagger_q_alpha a_p_alpha
                    if (spin_string_alpha_aaaa.create(q, sign_pq)) {

                        for (size_t r = 0; r < K; r++) {  // r loops over SOs
                            int sign_pqr = sign_pq;  // sign of the operator a_r_alpha a^dagger_q_alpha a_p_alpha
                            if (spin_string_alpha_aaaa.annihilate(r, sign_pqr)) {

                                for (size_t s = 0; s < K; s++) {  // s loops over SOs
                                    int sign_pqrs = sign_pqr;  // sign of the operator a^dagger_s_alpha a_r_alpha a^dagger_q_alpha a_p_alpha
                                    if (spin_string_alpha_aaaa.create(s, sign_pqrs)) {
                                        size_t J_alpha = fock_space_alpha.getAddress(spin_string_alpha_aaaa);  // the address of the string J_alpha that couples to I_alpha

                                        for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {  // I_beta loops over all beta addresses
                                            matvec(I_alpha*dim_beta + I_beta) += 0.5 * hamiltonian_parameters.get_g().get(p,q,r,s) * sign_pqrs * x(J_alpha*dim_beta + I_beta);
                                        }

                                        spin_string_alpha_aaaa.annihilate(s);  // undo the previous creation
                                    }
                                }  // loop over s

                                spin_string_alpha_aaaa.create(r);  // undo the previous annihilation
                            }
                        }  // loop over r

                        spin_string_alpha_aaaa.annihilate(q);  // undo the previous creation
                    }
                }  // loop over q

                spin_string_alpha_aaaa.create(p);  // undo the previous creation
            }
        }  // loop over p
    }  // loop over I_alpha


    // ALPHA-ALPHA-BETA-BETA (and BETA-BETA-ALPHA-ALPHA)
    ONV spin_string_alpha_aabb = fock_space_alpha.get_ONV(0);  // spin string with address 0
    for (size_t I_alpha = 0; I_alpha < dim_alpha; I_alpha++) {  // I_alpha loops over all addresses of alpha spin strings
        if (I_alpha > 0) {
            fock_space_alpha.setNext(spin_string_alpha_aabb);
        }

        for (size_t p = 0; p < K; p++) {  // p loops over SOs
            int sign_p = 1;  // sign of the operator a_p_alpha
            if (spin_string_alpha_aabb.annihilate(p, sign_p)) {

                for (size_t q = 0; q < K; q++) {
                    int sign_pq = sign_p;  // sign of the operator a^dagger_q_alpha a_p_alpha
                    if (spin_string_alpha_aabb.create(q, sign_pq)) {
                        size_t J_alpha = fock_space_alpha.getAddress(spin_string_alpha_aabb);  // the address of the spin string that couples to I_alpha

                        ONV spin_string_beta_aabb = fock_space_beta.get_ONV (0); // spin string with address 0
                        for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {  // I_beta loops over all addresses of beta spin strings
                            if (I_beta > 0) {
                                fock_space_beta.setNext(spin_string_beta_aabb);
                            }

                            for (size_t r = 0; r < K; r++) {  // r loops over SOs
                                int sign_r = 1;  // sign of the operator a_r_beta
                                if (spin_string_beta_aabb.annihilate(r, sign_r)) {

                                    for (size_t s = 0; s < K; s++) {  // s loops over SOs
                                        int sign_rs = sign_r;  // sign of the operato a^dagger_s_beta a_r_beta
                                        if (spin_string_beta_aabb.create(s, sign_rs)) {
                                            size_t J_beta = fock_space_beta.getAddress(spin_string_beta_aabb);  // the address of the spin string that couples to I_beta

                                            matvec(I_alpha*dim_beta + I_beta) += hamiltonian_parameters.get_g().get(p,q,r,s) * sign_pq * sign_rs * x(J_alpha*dim_beta + J_beta);  // alpha addresses are major

                                            spin_string_beta_aabb.annihilate(s);  // undo the previous creation
                                        }
                                    }  // loop over r

                                    spin_string_beta_aabb.create(r);  // undo the previous annihilation
                                }
                            }  // loop over r


                        }  // I_beta loop

                        spin_string_alpha_aabb.annihilate(q);  // undo the previous creation
                    }
                }  // loop over q

                spin_string_alpha_aabb.create(p);  // undo the previous annihilation
            }
        }  // loop over p
    }  // loop over I_alpha


    // BETA-BETA-BETA-BETA
    ONV spin_string_beta_bbbb = fock_space_beta.get_ONV(0);  // spin string with address 0
    for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {  // I_beta loops over all addresses of beta spin strings
        if (I_beta > 0) {
            fock_space_beta.setNext(spin_string_beta_bbbb);
        }

        for (size_t p = 0; p < K; p++) {  // p loops over SOs
            int sign_p = 1;  // sign of the operator a_p_beta
            if (spin_string_beta_bbbb.annihilate(p, sign_p)) {

                for (size_t q = 0; q < K; q++) {  // q loops over SOs
                    int sign_pq = sign_p;  // sign of the operator a^dagger_q_beta a_p_beta
                    if (spin_string_beta_bbbb.create(q, sign_pq)) {

                        for (size_t r = 0; r < K; r++) {  // r loops over SOs
                            int sign_pqr = sign_pq;  // sign of the operator a_r_beta a^dagger_q_beta a_p_beta
                            if (spin_string_beta_bbbb.annihilate(r, sign_pqr)) {

                                for (size_t s = 0; s < K; s++) {  // s loops over SOs
                                    int sign_pqrs = sign_pqr;  // sign of the operator a^dagger_s_beta a_r_beta a^dagger_q_beta a_p_beta
                                    if (spin_string_beta_bbbb.create(s, sign_pqrs)) {
                                        size_t J_beta = fock_space_beta.getAddress(spin_string_beta_bbbb);  // the address of the string J_beta that couples to I_beta

                                        for (size_t I_alpha = 0; I_alpha < dim_alpha; I_alpha++) {  // I_beta loops over all beta addresses
                                            matvec(I_alpha*dim_beta + I_beta) += 0.5 * hamiltonian_parameters.get_g().get(p,q,r,s) * sign_pqrs * x(I_alpha*dim_beta + J_beta);
                                        }

                                        spin_string_beta_bbbb.annihilate(s);  // undo the previous creation
                                    }
                                }  // loop over s

                                spin_string_beta_bbbb.create(r);  // undo the previous annihilation
                            }
                        }  // loop over r

                        spin_string_beta_bbbb.annihilate(q);  // undo the previous creation
                    }
                }  // loop over q

                spin_string_beta_bbbb.create(p);  // undo the previous creation
            }
        }  // loop over p
    }  // loop over I_beta

    return matvec;
}


/**
 *  @return the diagonal of the matrix representation of the Hamiltonian given @param hamiltonian_parameters
 */
Eigen::VectorXd FCI::calculateDiagonal(const HamiltonianParameters& hamiltonian_parameters) {

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
    // TODO: move this to libwint
    Eigen::MatrixXd k_SO = hamiltonian_parameters.get_h().get_matrix_representation();
    for (size_t p = 0; p < K; p++) {
        for (size_t q = 0; q < K; q++) {
            for (size_t r = 0; r < K; r++) {
                k_SO(p,q) -= 0.5 * hamiltonian_parameters.get_g().get(p,r,r,q);
            }
        }
    }

    ONV spin_string_alpha = fock_space_alpha.get_ONV(0);
    for (size_t Ia = 0; Ia < dim_alpha; Ia++) {  // Ia loops over addresses of alpha spin strings

        ONV spin_string_beta = fock_space_beta.get_ONV(0);
        for (size_t Ib = 0; Ib < dim_beta; Ib++) {  // Ib loops over addresses of beta spin strings

            for (size_t p = 0; p < K; p++) {  // p loops over SOs

                if (spin_string_alpha.isOccupied(p)) {  // p is in Ia
                    diagonal(Ia * dim_beta + Ib) += k_SO(p, p);

                    for (size_t q = 0; q < K; q++) {  // q loops over SOs
                        if (spin_string_alpha.isOccupied(q)) {  // q is in Ia
                            diagonal(Ia * dim_beta + Ib) += 0.5 * hamiltonian_parameters.get_g().get(p, p, q, q);
                        } else {  // q is not in I_alpha
                            diagonal(Ia * dim_beta + Ib) += 0.5 * hamiltonian_parameters.get_g().get(p, q, q, p);
                        }

                        if (spin_string_beta.isOccupied(q)) {  // q is in Ib
                            diagonal(Ia * dim_beta + Ib) += hamiltonian_parameters.get_g().get(p, p, q, q);
                        }
                    }  // q loop
                }


                if (spin_string_beta.isOccupied(p)) {  // p is in Ib
                    diagonal(Ia * dim_beta + Ib) += k_SO(p, p);


                    for (size_t q = 0; q < K; q++) {  // q loops over SOs
                        if (spin_string_beta.isOccupied(q)) {  // q is in Ib
                            diagonal(Ia * dim_beta + Ib) += 0.5 * hamiltonian_parameters.get_g().get(p, p, q, q);

                        } else {  // q is not in I_beta
                            diagonal(Ia * dim_beta + Ib) += 0.5 * hamiltonian_parameters.get_g().get(p, q, q, p);
                        }
                    }  // q loop
                }

            }  // p loop

            if (Ib < dim_beta - 1) {  // prevent last permutation to occur
                fock_space_beta.setNext(spin_string_beta);
            }
        }  // beta address (Ib) loop

        if (Ia < dim_alpha - 1) {  // prevent last permutation to occur
            fock_space_alpha.setNext(spin_string_alpha);
        }
    }  // alpha address (Ia) loop

    return diagonal;
}



}  // namespace GQCG
