#include "Drivers/FullConfigurationInteractionDriver.hpp"

namespace GQCP {

    FullConfigurationInteractionDriver::FullConfigurationInteractionDriver(std::string xyz_filename,
                                                                           std::string basis_set, size_t num_alpha,
                                                                           size_t num_beta) {

        molecule = std::make_shared<GQCP::Molecule>(GQCP::Molecule::Readxyz(xyz_filename));
        auto mol_ham_par = GQCP::HamiltonianParameters<double>::Molecular(*molecule, basis_set);  // in the AO basis
        mol_ham_par.LowdinOrthonormalize();  // now in the Löwdin basis

        // Solve the FCI eigenvalue problem using the dense algorithm
        auto K = mol_ham_par.get_K();
        GQCP::ProductFockSpace fock_space(K, num_alpha, num_beta);
        GQCP::FCI fci(fock_space);
        solver = std::make_shared<GQCP::CISolver>(fci, mol_ham_par);
        GQCP::DenseSolverOptions dense_solver_options;
        solver->solve(dense_solver_options);
    }

    double FullConfigurationInteractionDriver::get_energy() {
        double fci_energy = this->solver->get_eigenpair().get_eigenvalue();
        double internuclear_repulsion_energy = this->molecule->calculateInternuclearRepulsionEnergy();
        return fci_energy + internuclear_repulsion_energy;
    }
}