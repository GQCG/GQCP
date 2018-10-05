#ifndef GQCG_FOCKSPACE_HPP
#define GQCG_FOCKSPACE_HPP


#include "FockSpace/BaseFockSpace.hpp"
#include <boost/numeric/conversion/converter.hpp>
#include <boost/math/special_functions.hpp>


namespace GQCG {


class FockSpace: public GQCG::BaseFockSpace {
private:
    const size_t N;  // number of electrons
    const size_t dim;  // dimension of the Fock space
    Matrixu vertex_weights;


    // PRIVATE METHODS
    /**
     *  In-place permute the representation, giving the next bitstring permutation in reverse lexical ordering.
     *
     *      Examples:
     *          011 -> 101
     *          101 -> 110
     */
    size_t ulongNextPermutation(size_t representation) ;



public:
    // CONSTRUCTORS
    /**
     *  Constructor given a @param K, N
     */
    explicit FockSpace(size_t K, size_t N);


    // DESTRUCTOR
    ~FockSpace() override = default;

    // GETTERS
    /**
     *  @return weights as size_t from the vertex_weight matrix associated with the ONVs in the Fock space
     */
    size_t get_vertex_weights(size_t p, size_t m) const { return this->vertex_weights[p][m];}
    Matrixu get_vertex_weights() const { return this->vertex_weights;}

    size_t get_dimension(){ return dim;}

    /**
     *  @return the address (i.e. the ordering number) of the @param onv in reverse lexical ordering, in the fock space.
     */
    size_t get_address(ONV &onv);


    // STATIC PUBLIC METHODS
    /**
     *  Given a number of spatial orbitals @param K and a number of electrons  @param N, @return the dimension of
     *  the Fock space.
     */
    static size_t calculateDimension(size_t K, size_t N);


    // OVERRIDDEN PUBLIC METHODS
    /**
     *  @return ONV with the corresponding address in the considered space
     */
    ONV get_ONV(size_t address) override;

    /**
     *  sets @param ONV to the next ONV in the space
     *  performs the ulongNextPermutation() function
     *  and updates the corresponding occupation indexes
     */
    void setNext(ONV &onv) override;


    // FRIEND CLASSES
    friend class DOCI;
    friend class FCI;
};

typedef std::shared_ptr<GQCG::FockSpace> FockSpace_sptr;

}  // namespace GQCG


#endif //GQCG_FOCKSPACE_HPP
