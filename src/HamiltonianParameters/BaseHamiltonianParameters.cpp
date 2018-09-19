#include "HamiltonianParameters/BaseHamiltonianParameters.hpp"


namespace GQCG {


/*
 *  CONSTRUCTORS
 */

/**
 *  Constructor based on a given @param ao_basis
 */
BaseHamiltonianParameters::BaseHamiltonianParameters(const GQCG::AOBasis& ao_basis) :
    ao_basis (ao_basis)
{}



/*
 *  DESTRUCTOR
 */

/**
 *  Provide a pure virtual destructor to make the class abstract
 */
BaseHamiltonianParameters::~BaseHamiltonianParameters() {}



}  // namespace GQCG
