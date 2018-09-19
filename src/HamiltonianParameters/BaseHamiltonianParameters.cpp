#include "HamiltonianParameters/BaseHamiltonianParameters.hpp"


namespace GQCG {


/*
 *  CONSTRUCTORS
 */

/**
 *  Constructor based on a given @param ao_basis_ptr
 */
BaseHamiltonianParameters::BaseHamiltonianParameters(std::shared_ptr<GQCG::AOBasis> ao_basis_ptr) :
    ao_basis_ptr (std::move(ao_basis_ptr))
{}



/*
 *  DESTRUCTOR
 */

/**
 *  Provide a pure virtual destructor to make the class abstract
 */
BaseHamiltonianParameters::~BaseHamiltonianParameters() {}



}  // namespace GQCG
