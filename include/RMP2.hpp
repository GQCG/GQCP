#ifndef RMP2_hpp
#define RMP2_hpp

#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "Molecule.hpp"
#include "RHF/RHF.hpp"


namespace GQCG {


/**
 *  @return the RMP2 energy correction based on given @param Hamiltonian parameters ham_par, a given @param molecule and a converged solution @param rhf to the RHF SCF equations
 */
double calculateRMP2EnergyCorrection(const GQCG::HamiltonianParameters& ham_par, const GQCG::Molecule& molecule, const GQCG::RHF& rhf);


}  // namespace GQCG



#endif /* RMP2_hpp */
