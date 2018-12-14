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
#include "RDM/DOCIRDMBuilder.hpp"


namespace GQCP {


/*
 *  CONSTRUCTOR
 */
DOCIRDMBuilder::DOCIRDMBuilder(const FockSpace& fock_space) :
    fock_space (fock_space)
{}


/*
 *  OVERRIDDEN PUBLIC METHODS
 */

/**
 *  @param x        the coefficient vector representing the DOCI wave function
 *
 *  @return all 1-RDMs given a coefficient vector
 */
OneRDMs DOCIRDMBuilder::calculate1RDMs(const Eigen::VectorXd& x) const {
    // The formulas for the DOCI 1-RDMs can be found in (https://github.com/lelemmen/electronic_structure)

    size_t K = this->fock_space.get_K();
    size_t dim = this->fock_space.get_dimension();

    // For DOCI, one rdm covers both spins
    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(K, K);

    // Create the first ONV (with address 0). In DOCI, the Fock space for alpha and beta is equal so we just use one
    ONV onv = this->fock_space.makeONV(0);   


    for (size_t I = 0; I < dim; I++) {  // I loops over all the addresses of the spin strings

        for (size_t e1 = 0; e1 < this->fock_space.get_N(); e1++) {  // e1 (electron 1) loops over the (number of) electrons
            size_t p = onv.get_occupation_index(e1);  // retrieve the index of the orbital the electron occupies
            double c_I = x(I);  // coefficient of the I-th basis vector
            D(p,p) += 2*std::pow(c_I, 2);
        }
        
        if (I < dim-1) {
            this->fock_space.setNextONV(onv);
        }
    }

    OneRDM one_rdm (D);
    return OneRDMs (one_rdm);
}


/**
 *  @param x        the coefficient vector representing the DOCI wave function
 *
 *  @return all 2-RDMs given a coefficient vector
 */
TwoRDMs DOCIRDMBuilder::calculate2RDMs(const Eigen::VectorXd& x) const {

    // The formulas for the DOCI 2-RDMs can be found in (https://github.com/lelemmen/electronic_structure)

    size_t K = this->fock_space.get_K();
    size_t dim = this->fock_space.get_dimension();
    
    Eigen::Tensor<double, 4> d_aaaa (K, K, K, K);
    d_aaaa.setZero();
    Eigen::Tensor<double, 4> d_aabb (K, K, K, K);
    d_aabb.setZero();

    // Create the first ONV (with address 0). In DOCI, the Fock space for alpha and beta is equal so we just use one
    ONV onv = this->fock_space.makeONV(0);   


    for (size_t I = 0; I < dim; I++) {  // I loops over all the addresses of the spin strings
        for (size_t p = 0; p < K; p++) {  // p loops over SOs
            if (onv.annihilate(p)) {  // if p is occupied in I

                double c_I = x(I);  // coefficient of the I-th basis vector
                double c_I_2 = std::pow(c_I, 2);  // square of c_I

                d_aabb(p,p,p,p) += c_I_2;

                for (size_t q = 0; q < p; q++) {  // q loops over SOs with an index smaller than p
                    if (onv.create(q)) {  // if q is not occupied in I
                        size_t J = this->fock_space.getAddress(onv);  // the address of the coupling string
                        double c_J = x(J);  // coefficient of the J-th basis vector

                        d_aabb(p,q,p,q) += c_I * c_J;
                        d_aabb(q,p,q,p) += c_I * c_J;  // since we're looping for q < p

                        onv.annihilate(q);  // reset the spin string after previous creation on q
                    }

                    else {  // if q is occupied in I
                        d_aaaa(p,p,q,q) += c_I_2;
                        d_aaaa(q,q,p,p) += c_I_2;  // since we're looping for q < p

                        d_aaaa(p,q,q,p) -= c_I_2;
                        d_aaaa(q,p,p,q) -= c_I_2;  // since we're looping for q < p

                        d_aabb(p,p,q,q) += c_I_2;
                        d_aabb(q,q,p,p) += c_I_2;  // since we're looping for q < p
                    }
                }
                onv.create(p);  // reset the spin string after previous annihilation on p
            }
        }

        if ( I < dim-1) {
            this->fock_space.setNextONV(onv);
        }
    }

    // For DOCI, we have additional symmetries (two_rdm_aaaa = two_rdm_bbbb, two_rdm_aabb = two_rdm_bbaa)
    TwoRDM two_rdm_aaaa (d_aaaa);
    TwoRDM two_rdm_aabb (d_aabb);
    TwoRDM two_rdm_bbaa (d_aabb);
    TwoRDM two_rdm_bbbb (d_aaaa);
    return TwoRDMs (two_rdm_aaaa, two_rdm_aabb, two_rdm_bbaa, two_rdm_bbbb);
}


}  // namespace GQCP
