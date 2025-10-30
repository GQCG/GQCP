#include "Basis/SpinorBasis/RSpinOrbitalBasis.hpp"
#include "QCMethod/HF/RHF/DiagonalRHFFockMatrixObjective.hpp"
#include "QCMethod/HF/RHF/RHF.hpp"
#include "QCMethod/HF/RHF/RHFSCFSolver.hpp"
#include "QCModel/HF/RHF.hpp"

int main() {

    // Create an STO-3G basisset on (a toy geometry of) H2O.
    const GQCP::Nucleus h1 {1, 0.0, 0.0, 0.0};
    const GQCP::Nucleus o {8, 1.0, 0.0, 0.0};
    const GQCP::Nucleus h2 {1, 2.0, 0.0, 0.0};
    const GQCP::Molecule molecule {{h1, o, h2}};

    // Perform an RHF calculation.
    const GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "STO-3G"};
    auto hamiltonian = spin_orbital_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));  // In an AO basis.

    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(molecule.numberOfElectrons(), hamiltonian, spin_orbital_basis.overlap());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective {hamiltonian};

    const auto rhf_qc_structure = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment);
    const auto rhf_parameters = rhf_qc_structure.groundStateParameters();
    const auto rhf_energy = rhf_qc_structure.groundStateEnergy();

    // Determine the RHF energy through the expectation value of the Hamiltonian, and check the result.
    // Do the calculations in the RHF MO basis, in order to check the implementation of the RHF density matrices in MO basis.
    hamiltonian.transform(rhf_parameters.expansion());
    const auto D_MO = rhf_parameters.calculateOrthonormalBasis1DM();
    const auto d_MO = rhf_parameters.calculateOrthonormalBasis2DM();
    const double expectation_value = hamiltonian.calculateExpectationValue(D_MO, d_MO);

    // obtain F (established)
    const auto& D = rhf_environment.density_matrices.back();  // The most recent density matrix.
    const auto F = GQCP::QCModel::RHF<double>::calculateScalarBasisFockMatrix(D, rhf_environment.sq_hamiltonian);

    std::cout << "F" << std::endl;
    std::cout << F.parameters() << std::endl;

    // obtain J and K (new functions)
    const auto J = GQCP::QCModel::RHF<double>::calculateScalarBasisDirectMatrix(D, rhf_environment.sq_hamiltonian);
    const auto K = GQCP::QCModel::RHF<double>::calculateScalarBasisExchangeMatrix(D, rhf_environment.sq_hamiltonian);

    std::cout << "J" << std::endl;
    std::cout << J.parameters() << std::endl;

    std::cout << "K" << std::endl;
    std::cout << K.parameters() << std::endl;

    return 0;
}