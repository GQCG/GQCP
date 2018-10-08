#include <FockSpace/FockSpace.hpp>


namespace GQCG {


/*
 * PRIVATE METHODS
 */

/**
 *  In-place permute the unsigned representation of the @param ONV, giving the next bitstring permutation in reverse lexical ordering.
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
 *  Constructor given a @param K (spatial orbitals), N (electrons)
 *  on which the dimensions of the Fock space are based
 */

FockSpace::FockSpace(size_t K, size_t N) :
        BaseFockSpace(K), N(N), dim(FockSpace::calculateDimension(K,N)) {
    // Create a zero matrix of dimensions (K+1)x(N+1)
    this->vertex_weights = GQCG::Matrixu(this->K + 1, GQCG::Vectoru(this->N + 1, 0));

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
 *  Given a number of spatial orbitals @param K
 *  and a number of electrons  @param N,
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
 *  @return the ONV with the corresponding address in the considered space
 */
ONV FockSpace::get_ONV(size_t address) {
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
    return ONV(this->K, this->N, representation);
}


/**
 *  sets @param ONV to the next ONV in the space
 *  performs the ulongNextPermutation() function
 *  and updates the corresponding occupation indexes
 *  of the ONV occupation array
 */
void FockSpace::setNext(ONV& onv) {
    onv.set_representation(ulongNextPermutation(onv.unsigned_representation));
}


/**
 *  @return the Fock space address (i.e. the ordering number) of the @param onv in reverse lexical ordering, in the fock space.
 */
size_t FockSpace::getAddress(ONV& onv) {
    // An implementation of the formula in Helgaker, starting the addressing count from zero
    size_t address = 0;
    size_t electron_count = 0;  // counts the number of electrons in the spin string up to orbital p
    unsigned long unsigned_onv = onv.unsigned_representation;  // copy the unsigned_representation of the onv

    while(unsigned_onv != 0) {  // we will remove the least significant bit each loop, we are finished when no bits are left
        size_t p = __builtin_ctzl(unsigned_onv);  // p is the orbital index counter (starting from 1)
        electron_count++;  // each bit is an electron hence we add it up to the electron count
        address += get_vertex_weights(p , electron_count);
        unsigned_onv ^= unsigned_onv & -unsigned_onv;  // flip the least significant bit
    }
    return address;
}


}  // namespace GQCG
