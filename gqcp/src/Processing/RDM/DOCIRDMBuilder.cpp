// This file is part of GQCG-GQCP.
//
// Copyright (C) 2017-2020  the GQCG developers
//
// GQCG-GQCP is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GQCG-GQCP is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-GQCP.  If not, see <http://www.gnu.org/licenses/>.

#include "Processing/RDM/DOCIRDMBuilder.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param onv_basis            the seniority-zero ONV basis
 */
DOCIRDMBuilder::DOCIRDMBuilder(const SeniorityZeroONVBasis& onv_basis) :
    onv_basis {onv_basis} {}


/*
 *  OVERRIDDEN PUBLIC METHODS
 */

/**
 *  @param x        the coefficient vector representing the DOCI wave function
 *
 *  @return the 1-RDMs given a coefficient vector
 */
OneRDMs<double> DOCIRDMBuilder::calculate1RDMs(const VectorX<double>& x) const {

    // Prepare some variables.
    const auto K = this->onv_basis.numberOfSpatialOrbitals();
    const auto dimension = this->onv_basis.dimension();

    // For seniority-zero linear expansions, one DM covers both alpha and beta spins.
    OneRDM<double> D = OneRDM<double>::Zero(K, K);

    // Create the first ONV (with address 0). In DOCI, the ONV basis for alpha and beta is equal, so we can use the proxy ONV basis.
    const auto onv_basis_proxy = this->onv_basis.proxy();
    SpinUnresolvedONV onv = onv_basis_proxy.makeONV(0);
    for (size_t I = 0; I < dimension; I++) {  // I loops over all the addresses of the doubly-occupied ONVs

        for (size_t e1 = 0; e1 < onv_basis_proxy.get_N(); e1++) {  // e1 (electron 1) loops over the number of electrons
            const size_t p = onv.occupationIndexOf(e1);         // retrieve the index of the orbital the electron occupies
            const double c_I = x(I);                               // coefficient of the I-th basis vector

            D(p, p) += 2 * std::pow(c_I, 2);
        }

        if (I < dimension - 1) {  // prevent the last permutation from occurring
            onv_basis_proxy.setNextONV(onv);
        }
    }

    return OneRDMs<double>(D);
}


/**
 *  @param x        the coefficient vector representing the DOCI wave function
 *
 *  @return the 2-RDMs given a coefficient vector
 */
TwoRDMs<double> DOCIRDMBuilder::calculate2RDMs(const VectorX<double>& x) const {

    // Prepare some variables.
    const auto K = this->onv_basis.numberOfSpatialOrbitals();
    const auto dimension = this->onv_basis.dimension();


    // For seniority-zero linear expansions, we only have to calculate d_aaaa and d_aabb.
    TwoRDM<double> d_aaaa(K);
    d_aaaa.setZero();

    TwoRDM<double> d_aabb(K);
    d_aabb.setZero();


    // Create the first ONV (with address 0). In DOCI, the ONV basis for alpha and beta is equal, so we can use the proxy ONV basis.
    const auto onv_basis_proxy = this->onv_basis.proxy();
    SpinUnresolvedONV onv = onv_basis_proxy.makeONV(0);
    for (size_t I = 0; I < dimension; I++) {  // I loops over all the addresses of the spin strings
        for (size_t p = 0; p < K; p++) {      // p loops over SOs
            if (onv.annihilate(p)) {          // if p is occupied in I

                const double c_I = x(I);                // coefficient of the I-th basis vector
                const double c_I_2 = std::pow(c_I, 2);  // square of c_I

                d_aabb(p, p, p, p) += c_I_2;

                for (size_t q = 0; q < p; q++) {                           // q loops over SOs with an index smaller than p
                    if (onv.create(q)) {                                   // if q is not occupied in I
                        const size_t J = onv_basis_proxy.getAddress(onv);  // the address of the coupling string
                        const double c_J = x(J);                           // coefficient of the J-th basis vector

                        d_aabb(p, q, p, q) += c_I * c_J;
                        d_aabb(q, p, q, p) += c_I * c_J;  // since we're looping for q < p

                        onv.annihilate(q);  // reset the spin string after previous creation on q
                    }

                    else {  // if q is occupied in I
                        d_aaaa(p, p, q, q) += c_I_2;
                        d_aaaa(q, q, p, p) += c_I_2;  // since we're looping for q < p

                        d_aaaa(p, q, q, p) -= c_I_2;
                        d_aaaa(q, p, p, q) -= c_I_2;  // since we're looping for q < p

                        d_aabb(p, p, q, q) += c_I_2;
                        d_aabb(q, q, p, p) += c_I_2;  // since we're looping for q < p
                    }
                }
                onv.create(p);  // reset the spin string after previous annihilation on p
            }
        }

        if (I < dimension - 1) {  // prevent the last permutation from occurring
            onv_basis_proxy.setNextONV(onv);
        }
    }

    // For seniority-zero linear expansions, we have additional symmetries (two_rdm_aaaa = two_rdm_bbbb, two_rdm_aabb = two_rdm_bbaa)
    return TwoRDMs<double>(d_aaaa, d_aabb, d_aabb, d_aaaa);
}


}  // namespace GQCP
