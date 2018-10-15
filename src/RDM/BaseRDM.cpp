#include "RDM/BaseRDM.cpp"


namespace GQCG {



/*
 *  CONSTRUCTORS
 */

/**
 *  Constructor based on a given @param dimension
 */
BaseRDM::BaseRDM(size_t dimension) :
    dim (dimension)
{}


/*
 *  DESTRUCTOR
 */

/**
 *  Provide a pure virtual destructor to make the class abstract
 */
BaseRDM::~BaseRDM() {}


}  // namespace GQCG
