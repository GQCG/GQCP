#include "HamiltonianParameters/HamiltonianParameters_constructors.hpp"

#include "LibintCommunicator.hpp"


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
GQCG::HamiltonianParameters constructMolecularHamiltonianParameters(std::shared_ptr<GQCG::AOBasis> ao_basis_sptr) {

    // Calculate the integrals for the molecular Hamiltonian
    auto S = GQCG::LibintCommunicator::get().calculateOneElectronIntegrals(libint2::Operator::overlap, *ao_basis_sptr);
    auto T = GQCG::LibintCommunicator::get().calculateOneElectronIntegrals(libint2::Operator::kinetic, *ao_basis_sptr);
    auto V = GQCG::LibintCommunicator::get().calculateOneElectronIntegrals(libint2::Operator::nuclear, *ao_basis_sptr);
    auto H = T + V;
    
    auto g = GQCG::LibintCommunicator::get().calculateTwoElectronIntegrals(libint2::Operator::coulomb, *ao_basis_sptr);
    
    
    // Construct the initial transformation matrix: the identity matrix
    auto nbf = ao_basis_sptr->get_number_of_basis_functions();
    Eigen::MatrixXd C = Eigen::MatrixXd::Identity(nbf, nbf);
    
    
    return HamiltonianParameters(ao_basis_sptr, S, H, g, C);
}



}  // namespace GQCG
