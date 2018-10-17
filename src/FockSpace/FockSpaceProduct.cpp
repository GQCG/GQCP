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
        BaseFockSpace (K, FockSpaceProduct::calculateDimension(K, N_alpha, N_beta)),
        fock_space_alpha (FockSpace(K, N_alpha)),
        fock_space_beta (FockSpace(K, N_beta)),
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
size_t FockSpaceProduct::calculateDimension(size_t K, size_t N_alpha, size_t N_beta) {
    size_t alpha_dim = FockSpace::calculateDimension(K, N_alpha);
    size_t beta_dim = FockSpace::calculateDimension(K, N_beta);
    return boost::numeric::converter<double, size_t>::convert(beta_dim * alpha_dim);
}


}  // namespace GQCG
