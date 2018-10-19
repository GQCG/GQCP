#include "HamiltonianParameters/BaseHamiltonianParameters.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  Constructor based on a given @param ao_basis_sptr
 */
BaseHamiltonianParameters::BaseHamiltonianParameters(std::shared_ptr<GQCP::AOBasis> ao_basis) :
    ao_basis (std::move(ao_basis))
{}



/*
 *  DESTRUCTOR
 */

/**
 *  Provide a pure virtual destructor to make the class abstract
 */
BaseHamiltonianParameters::~BaseHamiltonianParameters() {}



}  // namespace GQCP
