#ifndef GQCG_BASERDM_HPP
#define GQCG_BASERDM_HPP


#include "common.hpp"


namespace GQCG {


/**
 *  A base class for the representation of reduced density matrices
 */
class BaseRDM {
protected:
    const size_t dim;  // dimension of the matrix representation of the operator


public:
    // CONSTRUCTORS
    /**
     *  Constructor based on a given @param dimension
     */
    explicit BaseRDM(size_t dimension);


    // DESTRUCTOR
    /**
     *  Provide a pure virtual destructor to make the class abstract
     */
    virtual ~BaseRDM() = 0;


    // GETTERS
    size_t get_dim() { return this->dim; }
};


}  // namespace GQCG


#endif  // GQCG_BASERDM_HPP
