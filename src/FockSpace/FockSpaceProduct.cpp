#include "FockSpace/FockSpaceProduct.hpp"


namespace GQCG {



/*
 *  CONSTRUCTORS
 */

/**
 *  Constructor given a @param K (spatial orbitals), N_alpha and N_beta (electrons)
 *  on which the dimensions of the Fock space are based
 */

FockSpaceProduct::FockSpaceProduct(size_t K, size_t N_alpha, size_t N_beta) :
        BaseFockSpace (K, FockSpace::calculateDimension(K, N_alpha, N_beta)),
        N_alpha (N_alpha),
        N_beta (N_beta)
{}


/*
 *  STATIC PUBLIC METHODS
 */

/**
 *  Given a number of spatial orbitals @param K
 *  and a number of alpha and beta electrons @param N_alpha, N_beta,
 *  @return the dimension of the Fock space
 */
size_t FockSpace::calculateDimension(size_t K, size_t N_alpha, size_t N_beta) {
    size_t alpha_dim = FockSpace::calculateDimension(size_t K, size_t N_alpha);
    size_t beta_dim = FockSpace::calculateDimension(size_t K, size_t N_beta)
    return boost::numeric::converter<double, size_t>::convert(beta_dim * alpha_dim);
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
 *  and updates the corresponding occupation indices
 *  of the ONV occupation vector
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
