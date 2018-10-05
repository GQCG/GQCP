#ifndef GQCG_BASEFOCKSPACE_HPP
#define GQCG_BASEFOCKSPACE_HPP


#include "ONV.hpp"
#include "common.hpp"

namespace GQCG {


class BaseFockSpace {
protected:
    const size_t K;  // number spatial orbitals


    // PROTECTED CONSTRUCTORS
    /**
     *  Protected constructor given a @param K
     */
    explicit BaseFockSpace(size_t K):K(K){};



public:
    // DESTRUCTOR
    /**
     *  Provide a pure virtual destructor to make the class abstract
     */
    virtual ~BaseFockSpace() = 0;


    // PURE VIRTUAL PUBLIC METHODS
    /**
     *  @return ONV with the corresponding address in the considered space
     */
    virtual ONV get_ONV(size_t address) = 0;

    /**
     *  sets @param ONV to the next ONV in the space
     */
    virtual void setNext(ONV &onv) = 0;


};


}  // namespace GQCG


#endif //GQCG_BASEFOCKSPACE_HPP
