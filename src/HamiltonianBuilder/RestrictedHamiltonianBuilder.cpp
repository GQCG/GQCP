#include <HamiltonianBuilder/RestrictedHamiltonianBuilder.hpp>


namespace GQCG {


/*
 *  PROTECTED CONSTRUCTORS
 */

/**
 *  Protected constructor given a @param hamiltonian_parameters
 */
RestrictedHamiltonianBuilder::RestrictedHamiltonianBuilder(HamiltonianParameters hamiltonian_parameters) :
    hamiltonian_parameters (hamiltonian_parameters)
{}



/*
 *  DESTRUCTOR
 */

/**
 *  Provide a pure virtual destructor to make the class abstract
 */
RestrictedHamiltonianBuilder::~RestrictedHamiltonianBuilder() {}



}  // namespace GQCG
