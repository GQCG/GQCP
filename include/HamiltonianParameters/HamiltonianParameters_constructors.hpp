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


/**
 *  @return HamiltonianParameters corresponding to the contents of an @param fcidump_file
 */
GQCG::HamiltonianParameters readFCIDUMPFile(const std::string& fcidump_file);
 
 


}  // namespace GQCG




#endif  // GQCG_HAMILTONIANPARAMETERS_CONSTRUCTORS_HPP
