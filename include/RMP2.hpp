#ifndef RMP2_hpp
#define RMP2_hpp

#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "Molecule.hpp"


namespace GQCG {


/**
 *  Calculate and @return the RMP2 energy correction based on given @param Hamiltonian parameters ham_par and a given @param molecule
 */
double calculateRMP2EnergyCorrection(const GQCG::HamiltonianParameters& ham_par, const GQCG::Molecule& molecule);


}  // namespace GQCG



#endif /* RMP2_hpp */
