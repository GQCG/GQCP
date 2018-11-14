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
#include <FockSpace/FockSpace.hpp>

#include "FockSpace/FockSpace.hpp"


namespace GQCP {


/*
 * PRIVATE METHODS
 */

/**
 *  @param representation       a representation of an ONV
 *
 *  @return the next bitstring permutation
 *
 *      Examples:
 *          011 -> 101
 *          101 -> 110
 */
size_t FockSpace::ulongNextPermutation(size_t representation) {

    // t gets this->representation's least significant 0 bits set to 1
    unsigned long t = representation | (representation - 1UL);

    // Next set to 1 the most significant bit to change,
    // set to 0 the least significant ones, and add the necessary 1 bits.
    return (t + 1UL) | (((~t & (t+1UL)) - 1UL) >> (__builtin_ctzl(representation) + 1UL));
}



/*
 *  CONSTRUCTORS
 */

/**
 *  @param K        the number of orbitals
 *  @param N        the number of electrons
 */
FockSpace::FockSpace(size_t K, size_t N) :
        BaseFockSpace(K, FockSpace::calculateDimension(K, N)),
        N (N)
{
    // Create a zero matrix of dimensions (K+1)x(N+1)
    this->vertex_weights = GQCP::Matrixu(this->K + 1, GQCP::Vectoru(this->N + 1, 0));

    // K=5   N=2
    // [ 0 0 0 ]
    // [ 0 0 0 ]
    // [ 0 0 0 ]
    // [ 0 0 0 ]
    // [ 0 0 0 ]
    // [ 0 0 0 ]


    // The largest (reverse lexical) string is the one that includes the first (K-N+1) vertices of the first column
    //      This is because every vertical move from (p,m) to (p+1,m+1) corresponds to "orbital p+1 is unoccupied".
    //      Therefore, the largest reverse lexical string is the one where the first (K-N) orbitals are unoccupied.
    //      This means that there should be (K-N) vertical moves from (0,0).
    // Therefore, we may only set the weights of first (K-N+1) vertices of the first column to 1.
    for (size_t p = 0; p < this->K - this->N + 1; p++) {
        this->vertex_weights[p][0] = 1;
    }

    // K=5   N=2
    // [ 1 0 0 ]
    // [ 1 0 0 ]
    // [ 1 0 0 ]
    // [ 1 0 0 ]
    // [ 0 0 0 ]
    // [ 0 0 0 ]


    // The recurrence relation for the vertex weights is as follows:
    //      Every element is the sum of the values of the element vertically above and the element left diagonally above.
    //      W(p,m) = W(p-1,m) + W(p-1,m-1)

    for (size_t m = 1; m < this->N + 1; m++) {
        for (size_t p = m; p < (this->K - this->N + m) + 1; p++) {
            this->vertex_weights[p][m] = this->vertex_weights[p - 1][m] + this->vertex_weights[p - 1][m - 1];
        }
    }

    // K=5   N=2
    // [ 1 0 0 ]
    // [ 1 1 0 ]
    // [ 1 2 1 ]
    // [ 1 3 3 ]
    // [ 0 4 6 ]
    // [ 0 0 10]
}



/*
 *  STATIC PUBLIC METHODS
 */

/**
 *  @param K        the number of orbitals
 *  @param N        the number of electrons
 *
 *  @return the dimension of the Fock space
 */
size_t FockSpace::calculateDimension(size_t K, size_t N) {
    auto dim_double = boost::math::binomial_coefficient<double>(static_cast<unsigned>(K), static_cast<unsigned>(N));
    return boost::numeric::converter<double, size_t>::convert(dim_double);
}



/*
 *  PUBLIC METHODS
 */

/**
 *  @param address      the address (i.e. the ordening number) of the ONV
 *
 *  @return the ONV with the corresponding address
 */
ONV FockSpace::get_ONV(size_t address) {
    ONV onv (this->K, this->N);
    this->set(onv, address);
    return onv;
}


/**
 *  Set the current ONV to the next ONV: performs ulongNextPermutation() and updates the corresponding occupation indices of the ONV occupation array
 *
 *  @param onv      the current ONV
 */
void FockSpace::setNext(ONV& onv) {
    onv.set_representation(ulongNextPermutation(onv.get_unsigned_representation()));
}


/**
 *  @param onv      the ONV
 *
 *  @return the address (i.e. the ordering number) of the given ONV
 */
size_t FockSpace::getAddress(const ONV& onv) {
    // An implementation of the formula in Helgaker, starting the addressing count from zero
    size_t address = 0;
    size_t electron_count = 0;  // counts the number of electrons in the spin string up to orbital p
    unsigned long unsigned_onv = onv.get_unsigned_representation();  // copy the unsigned_representation of the onv

    while(unsigned_onv != 0) {  // we will remove the least significant bit each loop, we are finished when no bits are left
        size_t p = __builtin_ctzl(unsigned_onv);  // p is the orbital index counter (starting from 1)
        electron_count++;  // each bit is an electron hence we add it up to the electron count
        address += get_vertex_weights(p , electron_count);
        unsigned_onv ^= unsigned_onv & -unsigned_onv;  // flip the least significant bit
    }
    return address;
}


/**
 *  Transform an ONV to one with corresponding to the given address
 *
 *  @param onv          the ONV
 *  @param address      the address to which the ONV will be set
 */
void FockSpace::set(ONV& onv, size_t address) const {

    size_t representation;
    if (this->N == 0) {
        representation = 0;
    }

    else {
        representation = 0;
        size_t m = this->N;  // counts the number of electrons in the spin string up to orbital p

        for (size_t p = this->K; p > 0; p--) {  // p is an orbital index
            size_t weight = get_vertex_weights(p-1, m);

            if (weight <= address) {  // the algorithm can move diagonally, so we found an occupied orbital
                address -= weight;
                representation |= ((1) << (p - 1));  // set the (p-1)th bit: see (https://stackoverflow.com/a/47990)

                m--;  // since we found an occupied orbital, we have one electron less
                if (m == 0) {
                    break;
                }
            }
        }
    }
    onv.set_representation(representation);
}

  
 *  Find the next unoccupied orbital in a given ONV,
 *  update the electron count and orbital index,
 *  and calculate a shift in address
 *  resulting from a difference between the initial vertex weights for the encountered occupied orbitals
 *  and the corrected vertex weights accounting for previously annihilated electrons
 *
 *  @param onv       the ONV for which we search the next unnocupied orbital
 *  @param q         the orbital index
 *  @param e         the electron count
 *  @param a         the annihilation count
 *
 *  @return the shift in address resulting from the difference in the corrected electron weights
 */

size_t FockSpace::shiftUntilNextUnoccupiedOrbital(const ONV& onv, size_t& q, size_t& e, size_t a) const {

    size_t address_shift = 0;
    // Test whether the current orbital index is occupied
    while (e < this->N && q == onv.get_occupied_index(e)) {

        // Take the difference of vertex weights for the encountered electron weights to that of a vertex weight path with "a" fewer electrons
        // +1 is added to the electron index, because of how the addressing scheme is arrayed.
        address_shift += this->get_vertex_weights(q, e + 1 - a) - this->get_vertex_weights(q, e + 1);

        // move to the next electron and orbital
        e++;
        q++;
    }

    return address_shift;
}



}  // namespace GQCP
