// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2019  the GQCG developers
// 
// GQCG-gqcp is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// GQCG-gqcp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-gqcp.  If not, see <http://www.gnu.org/licenses/>.
// 
#define BOOST_TEST_MODULE "atomic_decomposition"

#include <boost/test/unit_test.hpp>

#include "QCMethod/Applications/AtomicDecompositionParameters.hpp"

#include "Basis/SpinorBasis/RSpinorBasis.hpp"
#include "Processing/Properties/expectation_values.hpp"
#include "Processing/RDM/RDMCalculator.hpp"
#include "QCMethod/CI/CISolver.hpp"
#include "QCMethod/CI/HamiltonianBuilder/FCI.hpp"
#include "QCMethod/HF/DiagonalRHFFockMatrixObjective.hpp"
#include "QCMethod/HF/RHF.hpp"
#include "QCMethod/HF/RHFSCFSolver.hpp"
#include "Utilities/units.hpp"


BOOST_AUTO_TEST_CASE ( decomposition_BeH_cation_STO_3G_Nuclear ) {

    // Create the molecular Hamiltonian in an AO basis
    GQCP::Nucleus Be(4, 0.0, 0.0, 0.0);
    GQCP::Nucleus H(1, 0.0, 0.0, GQCP::units::angstrom_to_bohr(1.134));  // from CCCBDB, STO-3G geometry
    std::vector<GQCP::Nucleus> nuclei {Be, H};
    GQCP::Molecule BeH (nuclei, +1);
    GQCP::AtomicDecompositionParameters adp = GQCP::AtomicDecompositionParameters::Nuclear(BeH, "STO-3G");
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (BeH, "STO-3G");

    auto sq_hamiltonian = adp.get_molecular_hamiltonian_parameters();
    auto K = sq_hamiltonian.dimension();
    double repulsion = GQCP::Operator::NuclearRepulsion(BeH).value();

    // Create a plain RHF SCF solver and solve the SCF equations
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(BeH.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective (sq_hamiltonian);
    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters();

    const auto& T = rhf_parameters.coefficientMatrix();

    // Transform the sq_hamiltonian
    sq_hamiltonian.transform(T);

    // Create the FCI module
    GQCP::ProductONVBasis fock_space (K, BeH.numberOfElectrons() / 2, BeH.numberOfElectrons() / 2);  // dim = 441
    GQCP::FCI fci (fock_space);
    GQCP::CISolver ci_solver (fci, sq_hamiltonian);

    // Solve Davidson
    GQCP::DavidsonSolverOptions solver_options (fock_space.HartreeFockExpansion());
    ci_solver.solve(solver_options);

    // Retrieve the eigenpair
    auto fci_coeff = ci_solver.get_eigenpair().get_eigenvector();
    auto fci_energy = ci_solver.get_eigenpair().get_eigenvalue();

    // Decomposition calculations require RDM-tracing
    GQCP::RDMCalculator rdm_calc(fock_space);
    rdm_calc.set_coefficients(fci_coeff);

    auto one_rdm = rdm_calc.calculate1RDMs().one_rdm;
    auto two_rdm = rdm_calc.calculate2RDMs().two_rdm;

    // Transform rdms to the AObasis
    GQCP::OneRDM<double> ao_one_rdm = one_rdm;
    GQCP::TwoRDM<double> ao_two_rdm = two_rdm;

    ao_one_rdm.basisTransformInPlace(T.adjoint());
    ao_two_rdm.basisTransformInPlace(T.adjoint());

    double self_energy_a = adp.get_net_atomic_parameters()[0].calculateExpectationValue(ao_one_rdm, ao_two_rdm);
    double self_energy_b = adp.get_net_atomic_parameters()[1].calculateExpectationValue(ao_one_rdm, ao_two_rdm);
    double interaction_energy_ab = adp.get_interaction_parameters()[0].calculateExpectationValue(ao_one_rdm, ao_two_rdm) + repulsion;
    double total_energy_a = adp.get_atomic_parameters()[0].calculateExpectationValue(ao_one_rdm, ao_two_rdm) + repulsion / 2;
    double total_energy_b = adp.get_atomic_parameters()[1].calculateExpectationValue(ao_one_rdm, ao_two_rdm) + repulsion / 2;

    BOOST_CHECK(std::abs(total_energy_a + total_energy_b - fci_energy - repulsion) < 1.0e-10);
    BOOST_CHECK(std::abs(self_energy_a + self_energy_b + interaction_energy_ab - fci_energy - repulsion) < 1.0e-10);
    BOOST_CHECK(std::abs(self_energy_a + interaction_energy_ab / 2 - total_energy_a) < 1.0e-10);
    BOOST_CHECK(std::abs(self_energy_b + interaction_energy_ab / 2 - total_energy_b) < 1.0e-10);
}
