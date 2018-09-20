#ifndef GQCG_HAMILTONIANPARAMETERS_CONSTRUCTORS_HPP
#define GQCG_HAMILTONIANPARAMETERS_CONSTRUCTORS_HPP



#include <memory>

#include "AOBasis.hpp"
#include "HamiltonianParameters.hpp"



namespace GQCG {


/**
 *  @return HamiltonianParameters corresponding to the molecular Hamiltonian for the given @param ao_basis_sptr
 *
 *  The molecular Hamiltonian has
 *      - one-electron contributions:
 *          - kinetic
 *          - nuclear attraction
 *      - two-electron contributions:
 *          - Coulomb repulsion
 */
GQCG::HamiltonianParameters constructMolecularHamiltonianParameters(std::shared_ptr<GQCG::AOBasis> ao_basis_sptr);




}  // namespace GQCG




#endif  // GQCG_HAMILTONIANPARAMETERS_CONSTRUCTORS_HPP
