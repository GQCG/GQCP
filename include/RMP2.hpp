#ifndef RMP2_hpp
#define RMP2_hpp

#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "Molecule.hpp"
#include "RHF/RHF.hpp"


namespace GQCP {


/**
 *  @return the RMP2 energy correction based on given @param Hamiltonian parameters ham_par, a given @param molecule and a converged solution @param rhf to the RHF SCF equations
 */
double calculateRMP2EnergyCorrection(const GQCP::HamiltonianParameters& ham_par, const GQCP::Molecule& molecule, const GQCP::RHF& rhf);


}  // namespace GQCP



#endif /* RMP2_hpp */
