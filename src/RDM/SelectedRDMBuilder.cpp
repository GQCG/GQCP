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
#include "RDM/SelectedRDMBuilder.hpp"

#include "FockSpace/Configuration.hpp"


namespace GQCP {


/*
 *  CONSTRUCTOR
 */
SelectedRDMBuilder::SelectedRDMBuilder(const SelectedFockSpace& fock_space) :
    fock_space (fock_space)
{}


/*
 *  OVERRIDDEN PUBLIC METHODS
 */

/**
 *  @param x        the coefficient vector representing the 'selected' wave function
 *
 *  @return all 1-RDMs given a coefficient vector
 */
OneRDMs SelectedRDMBuilder::calculate1RDMs(const Eigen::VectorXd& x) const {

    // Initialize as zero matrices
    size_t K = this->fock_space.get_K();

    Eigen::MatrixXd D_aa = Eigen::MatrixXd::Zero(K, K);
    Eigen::MatrixXd D_bb = Eigen::MatrixXd::Zero(K, K);
    
    
    size_t dim = fock_space.get_dimension();

    for (size_t I = 0; I < dim; I++) {  // loop over all addresses (1)
        Configuration configuration_I = this->fock_space.get_configuration(I);
        ONV alpha_I = configuration_I.onv_alpha;
        ONV beta_I = configuration_I.onv_beta;
        
        double c_I = x(I);


        // Calculate the diagonal of the 1-RDMs
        for (size_t p = 0; p < K; p++) {

            if (alpha_I.isOccupied(p)) {
                D_aa(p,p) += std::pow(c_I, 2);
            }

            if (beta_I.isOccupied(p)) {
                D_bb(p,p) += std::pow(c_I, 2);
            }
        }


        // Calculate the off-diagonal elements, by going over all other ONVs
        for (size_t J = I+1; J < dim; J++) {

            Configuration configuration_J = this->fock_space.get_configuration(J);
            ONV alpha_J = configuration_J.onv_alpha;
            ONV beta_J = configuration_J.onv_beta;

            double c_J = x(J);


            // 1 electron excitation in alpha (i.e. 2 differences), 0 in beta
            if ((alpha_I.countNumberOfDifferences(alpha_J) == 2) && (beta_I.countNumberOfDifferences(beta_J) == 0)) {

                // Find the orbitals that are occupied in one string, and aren't in the other
                size_t p = alpha_I.findDifferentOccupations(alpha_J)[0];  // we're sure that there is only 1 element in the std::vector<size_t>
                size_t q = alpha_J.findDifferentOccupations(alpha_I)[0];  // we're sure that there is only 1 element in the std::vector<size_t>

                // Calculate the total sign, and include it in the RDM contribution
                int sign = alpha_I.operatorPhaseFactor(p) * alpha_J.operatorPhaseFactor(q);
                D_aa(p,q) += sign * c_I * c_J;
                D_aa(q,p) += sign * c_I * c_J;
            }


            // 1 electron excitation in beta, 0 in alpha
            if ((alpha_I.countNumberOfDifferences(alpha_J) == 0) && (beta_I.countNumberOfDifferences(beta_J) == 2)) {

                // Find the orbitals that are occupied in one string, and aren't in the other
                size_t p = beta_I.findDifferentOccupations(beta_J)[0];  // we're sure that there is only 1 element in the std::vector<size_t>
                size_t q = beta_J.findDifferentOccupations(beta_I)[0];  // we're sure that there is only 1 element in the std::vector<size_t>

                // Calculate the total sign, and include it in the RDM contribution
                int sign = beta_I.operatorPhaseFactor(p) * beta_J.operatorPhaseFactor(q);
                D_bb(p,q) += sign * c_I * c_J;
                D_bb(q,p) += sign * c_I * c_J;
            }

        }  // loop over addresses J > I
    }  // loop over addresses I

    OneRDM one_rdm_aa (D_aa);
    OneRDM one_rdm_bb (D_bb);
    // The total 1-RDM is the sum of the spin components
    return OneRDMs (one_rdm_aa, one_rdm_bb);
}


/**
 *  @param x        the coefficient vector representing the 'selected' wave function
 *
 *  @return all 2-RDMs given a coefficient vector
 */
TwoRDMs SelectedRDMBuilder::calculate2RDMs(const Eigen::VectorXd& x) const {

    
    // Initialize as zero matrices
    size_t K = this->fock_space.get_K();

    size_t dim = fock_space.get_dimension();
    
    Eigen::Tensor<double, 4> d_aaaa (K,K,K,K);
    d_aaaa.setZero();
    Eigen::Tensor<double, 4> d_aabb (K,K,K,K);
    d_aabb.setZero();
    Eigen::Tensor<double, 4> d_bbaa (K,K,K,K);
    d_bbaa.setZero();
    Eigen::Tensor<double, 4> d_bbbb (K,K,K,K);
    d_bbbb.setZero();

    for (size_t I = 0; I < dim; I++) {  // loop over all addresses I

        Configuration configuration_I = this->fock_space.get_configuration(I);
        ONV alpha_I = configuration_I.onv_alpha;
        ONV beta_I = configuration_I.onv_beta;

        double c_I = x(I);
        
        for (size_t p = 0; p < K; p++) {

            // 'Diagonal' elements of the 2-RDM: aaaa and aabb
            if (alpha_I.isOccupied(p)) {
                for (size_t q = 0; q < K; q++) {
                    if (beta_I.isOccupied(q)) {
                        d_aabb(p,p,q,q) += std::pow(c_I, 2);
                    }

                    if (p != q) {  // can't create/annihilate the same orbital twice
                        if (alpha_I.isOccupied(q)) {
                            d_aaaa(p,p,q,q) += std::pow(c_I, 2);
                            d_aaaa(p,q,q,p) -= std::pow(c_I, 2);
                        }
                    }

                }  // loop over q
            }

            // 'Diagonal' elements of the 2-RDM: bbbb and bbaa
            if (beta_I.isOccupied(p)) {
                for (size_t q = 0; q < K; q++) {
                    if (alpha_I.isOccupied(q)) {
                        d_bbaa(p,p,q,q) += std::pow(c_I, 2);
                    }

                    if (p != q) {  // can't create/annihilate the same orbital twice
                        if (beta_I.isOccupied(q)) {
                            d_bbbb(p,p,q,q) += std::pow(c_I, 2);
                            d_bbbb(p,q,q,p) -= std::pow(c_I, 2);
                        }
                    }
                }  // loop over q
            }
        }  // loop over q


        for (size_t J = I+1; J < dim; J++) {

            Configuration configuration_J = this->fock_space.get_configuration(J);
            ONV alpha_J = configuration_J.onv_alpha;
            ONV beta_J = configuration_J.onv_beta;

            double c_J = x(J);

            // 1 electron excitation in alpha, 0 in beta
            if ((alpha_I.countNumberOfDifferences(alpha_J) == 2) && (beta_I.countNumberOfDifferences(beta_J) == 0)) {

                // Find the orbitals that are occupied in one string, and aren't in the other
                size_t p = alpha_I.findDifferentOccupations(alpha_J)[0];  // we're sure that there is only 1 element in the std::vector<size_t>
                size_t q = alpha_J.findDifferentOccupations(alpha_I)[0];  // we're sure that there is only 1 element in the std::vector<size_t>

                // Calculate the total sign
                int sign = alpha_I.operatorPhaseFactor(p) * alpha_J.operatorPhaseFactor(q);


                for (size_t r = 0; r < K; r++) {  // r loops over spatial orbitals

                    if (alpha_I.isOccupied(r) && alpha_J.isOccupied(r)) {  // r must be occupied on the left and on the right
                        if ((p != r) && (q != r)) {  // can't create or annihilate the same orbital
                            // Fill in the 2-RDM contributions
                            d_aaaa(p,q,r,r) += sign * c_I * c_J;
                            d_aaaa(r,q,p,r) -= sign * c_I * c_J;
                            d_aaaa(p,r,r,q) -= sign * c_I * c_J;
                            d_aaaa(r,r,p,q) += sign * c_I * c_J;

                            d_aaaa(q,p,r,r) += sign * c_I * c_J;
                            d_aaaa(q,r,r,p) -= sign * c_I * c_J;
                            d_aaaa(r,p,q,r) -= sign * c_I * c_J;
                            d_aaaa(r,r,q,p) += sign * c_I * c_J;
                        }
                    }

                    if (beta_I.isOccupied(r)) {  // beta_I == beta_J from the previous if-branch

                        // Fill in the 2-RDM contributions
                        d_aabb(p,q,r,r) += sign * c_I * c_J;
                        d_aabb(q,p,r,r) += sign * c_I * c_J;

                        d_bbaa(r,r,p,q) += sign * c_I * c_J;
                        d_bbaa(r,r,q,p) += sign * c_I * c_J;
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


                for (size_t r = 0; r < K; r++) {  // r loops over spatial orbitals

                    if (beta_I.isOccupied(r) && beta_J.isOccupied(r)) {  // r must be occupied on the left and on the right
                        if ((p != r) && (q != r)) {  // can't create or annihilate the same orbital
                            // Fill in the 2-RDM contributions
                            d_bbbb(p,q,r,r) += sign * c_I * c_J;
                            d_bbbb(r,q,p,r) -= sign * c_I * c_J;
                            d_bbbb(p,r,r,q) -= sign * c_I * c_J;
                            d_bbbb(r,r,p,q) += sign * c_I * c_J;

                            d_bbbb(q,p,r,r) += sign * c_I * c_J;
                            d_bbbb(q,r,r,p) -= sign * c_I * c_J;
                            d_bbbb(r,p,q,r) -= sign * c_I * c_J;
                            d_bbbb(r,r,q,p) += sign * c_I * c_J;
                        }
                    }

                    if (alpha_I.isOccupied(r)) {  // alpha_I == alpha_J from the previous if-branch

                        // Fill in the 2-RDM contributions
                        d_bbaa(p,q,r,r) += sign * c_I * c_J;
                        d_bbaa(q,p,r,r) += sign * c_I * c_J;

                        d_aabb(r,r,p,q) += sign * c_I * c_J;
                        d_aabb(r,r,q,p) += sign * c_I * c_J;
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

                // Calculate the total sign, and include it in the 2-RDM contribution
                int sign = alpha_I.operatorPhaseFactor(p) * alpha_J.operatorPhaseFactor(q) * beta_I.operatorPhaseFactor(r) * beta_J.operatorPhaseFactor(s);
                d_aabb(p,q,r,s) += sign * c_I * c_J;
                d_aabb(q,p,s,r) += sign * c_I * c_J;

                d_bbaa(r,s,p,q) += sign * c_I * c_J;
                d_bbaa(s,r,q,p) += sign * c_I * c_J;
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


                // Calculate the total sign, and include it in the 2-RDM contribution
                int sign = alpha_I.operatorPhaseFactor(p) * alpha_I.operatorPhaseFactor(r) * alpha_J.operatorPhaseFactor(q) * alpha_J.operatorPhaseFactor(s);
                d_aaaa(p,q,r,s) += sign * c_I * c_J;
                d_aaaa(p,s,r,q) -= sign * c_I * c_J;
                d_aaaa(r,q,p,s) -= sign * c_I * c_J;
                d_aaaa(r,s,p,q) += sign * c_I * c_J;

                d_aaaa(q,p,s,r) += sign * c_I * c_J;
                d_aaaa(s,p,q,r) -= sign * c_I * c_J;
                d_aaaa(q,r,s,p) -= sign * c_I * c_J;
                d_aaaa(s,r,q,p) += sign * c_I * c_J;
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


                // Calculate the total sign, and include it in the 2-RDM contribution
                int sign = beta_I.operatorPhaseFactor(p) * beta_I.operatorPhaseFactor(r) * beta_J.operatorPhaseFactor(q) * beta_J.operatorPhaseFactor(s);
                d_bbbb(p,q,r,s) += sign * c_I * c_J;
                d_bbbb(p,s,r,q) -= sign * c_I * c_J;
                d_bbbb(r,q,p,s) -= sign * c_I * c_J;
                d_bbbb(r,s,p,q) += sign * c_I * c_J;

                d_bbbb(q,p,s,r) += sign * c_I * c_J;
                d_bbbb(s,p,q,r) -= sign * c_I * c_J;
                d_bbbb(q,r,s,p) -= sign * c_I * c_J;
                d_bbbb(s,r,q,p) += sign * c_I * c_J;
            }

        }  // loop over all addresses J > I

    }  // loop over all addresses I

    TwoRDM two_rdm_aaaa (d_aaaa);
    TwoRDM two_rdm_aabb (d_aabb);
    TwoRDM two_rdm_bbaa (d_bbaa);
    TwoRDM two_rdm_bbbb (d_bbbb);
    return TwoRDMs (two_rdm_aaaa, two_rdm_aabb, two_rdm_bbaa, two_rdm_bbbb);
}


/**
 *  @param bra_indices      the indices of the orbitals that should be annihilated on the left (on the bra)
 *  @param ket_indices      the indices of the orbitals that should be annihilated on the right (on the ket)
 *  @param x                the coefficient vector representing the 'selected" wave function
 *
 *  @return an element of the N-RDM, as specified by the given bra and ket indices
 *
 *      calculateElement({0, 1}, {2, 1}) would calculate d^{(2)} (0, 1, 1, 2): the operator string would be a^\dagger_0 a^\dagger_1 a_2 a_1
 */
double SelectedRDMBuilder::calculateElement(const std::vector<size_t>& bra_indices, const std::vector<size_t>& ket_indices, const Eigen::VectorXd& x) const {
    throw std::runtime_error ("calculateElement is not implemented for SelectedRDMs");
}


}  // namespace GQCP
