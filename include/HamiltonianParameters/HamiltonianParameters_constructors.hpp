#ifndef GQCP_HAMILTONIANPARAMETERS_CONSTRUCTORS_HPP
#define GQCP_HAMILTONIANPARAMETERS_CONSTRUCTORS_HPP



#include <memory>

#include "AOBasis.hpp"
#include "HamiltonianParameters.hpp"



namespace GQCP {


/**
 *  @return HamiltonianParameters corresponding to the molecular Hamiltonian for the given @param ao_basis
 *
 *  The molecular Hamiltonian has
 *      - one-electron contributions:
 *          - kinetic
 *          - nuclear attraction
 *      - two-electron contributions:
 *          - Coulomb repulsion
 */
GQCP::HamiltonianParameters constructMolecularHamiltonianParameters(std::shared_ptr<GQCP::AOBasis> ao_basis);


/**
 *  @return HamiltonianParameters corresponding to the contents of an @param fcidump_file
 */
GQCP::HamiltonianParameters readFCIDUMPFile(const std::string& fcidump_file);
 

}  // namespace GQCP


#endif  // GQCP_HAMILTONIANPARAMETERS_CONSTRUCTORS_HPP
