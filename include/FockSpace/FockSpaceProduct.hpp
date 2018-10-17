#ifndef GQCG_FOCKSPACEPRODUCT_HPP
#define GQCG_FOCKSPACEPRODUCT_HPP


#include "FockSpace/BaseFockSpace.hpp"
#include "FockSpace/FockSpace.hpp"

#include <boost/numeric/conversion/converter.hpp>
#include <boost/math/special_functions.hpp>


namespace GQCG {


/**
 *  The product of two Fock spaces for a given set of orbitals and number of alpha and beta electrons.
 */
class FockSpaceProduct: public GQCG::BaseFockSpace {
private:
    const size_t N_alpha;  // number of alpha electrons
    const size_t N_beta;  // number of beta electrons

    FockSpace fock_space_alpha;
    FockSpace fock_space_beta;


public:
    // CONSTRUCTORS
    /**
     *  Constructor given a @param K (spatial orbitals), N_alpha and N_beta (electrons)
     *  on which the dimensions of the Fock space are based
     */
    FockSpaceProduct(size_t K, size_t N_alpha, size_t N_beta);


    // DESTRUCTORS
    ~FockSpaceProduct() override = default;


    // GETTERS
    size_t get_N_alpha() const { return this->N_alpha; }
    size_t get_N_beta() const { return this->N_beta; }
    FockSpace get_fock_space_alpha() const { return this->fock_space_alpha; }
    FockSpace get_fock_space_beta() const { return this->fock_space_beta; }


    // STATIC PUBLIC METHODS
    /**
     *  Given a number of spatial orbitals @param K
     *  and a number of alpha and beta electrons @param N_alpha, N_beta,
     *  @return the dimension of the Fock space
     */
    static size_t calculateDimension(size_t K, size_t N_alpha, size_t N_beta);
};


}  // namespace GQCG


#endif  // GQCG_FOCKSPACE_HPP
