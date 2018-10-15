#ifndef GQCG_FOCKSPACE_HPP
#define GQCG_FOCKSPACE_HPP


#include "FockSpace/BaseFockSpace.hpp"

#include <boost/numeric/conversion/converter.hpp>
#include <boost/math/special_functions.hpp>


namespace GQCG {


/**
 *  The full Fock space for a given set of orbitals and number of electrons
 *  where the ONVs and addresses are linked
 *  through a hashing function calculated with an addressing scheme.
 *  Implementation of the addressing scheme from :
 *      Molecular Electronic-Structure Theory (August 2000) by Trygve Helgaker, Poul Jorgensen, and Jeppe Olsen
 */
class FockSpace: public GQCG::BaseFockSpace {
private:
    const size_t N;  // number of electrons
    Matrixu vertex_weights;  // vertex_weights of the addressing scheme


    // PRIVATE METHODS
    /**
     *  @returns a permutation of the representation, giving the next bitstring permutation in reverse lexical ordering.
     *
     *      Examples:
     *          011 -> 101
     *          101 -> 110
     */
    size_t ulongNextPermutation(size_t representation);


public:
    // CONSTRUCTORS
    /**
     *  Constructor given a @param K (spatial orbitals), N (electrons)
     *  on which the dimensions of the Fock space are based
     */
    FockSpace(size_t K, size_t N);


    // DESTRUCTORS
    ~FockSpace() override = default;


    // GETTERS
    size_t get_vertex_weights(size_t p, size_t m) const { return this->vertex_weights[p][m]; }
    Matrixu get_vertex_weights() const { return this->vertex_weights; }


    // STATIC PUBLIC METHODS
    /**
     *  Given a number of spatial orbitals @param K
     *  and a number of electrons  @param N,
     *  @return the dimension of the Fock space
     */
    static size_t calculateDimension(size_t K, size_t N);


    // OVERRIDDEN PUBLIC METHODS
    /**
     *  @return the ONV with the corresponding address in the considered space
     */
    ONV get_ONV(size_t address) override;

    /**
     *  sets @param ONV to the next ONV in the space
     *  performs the ulongNextPermutation() function
     *  and updates the corresponding occupation indexes
     *  of the ONV occupation array
     */
    void setNext(ONV& onv) override;

    /**
     *  @return the Fock space address (i.e. the ordering number) of the @param onv in reverse lexical ordering, in the fock space.
     */
    size_t getAddress(ONV& onv) override;


    // FRIEND CLASSES
    friend class DOCI;
    friend class FCI;
};


}  // namespace GQCG


#endif  // GQCG_FOCKSPACE_HPP
