// This file is part of GQCG-GQCP.
//
// Copyright (C) 2017-2020  the GQCG developers
//
// GQCG-GQCP is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GQCG-GQCP is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-GQCP.  If not, see <http://www.gnu.org/licenses/>.

#define BOOST_TEST_MODULE "atomic_decomposition"

#include <boost/test/unit_test.hpp>

#include "Basis/SpinorBasis/RSpinorBasis.hpp"
#include "Mathematical/Optimization/Eigenproblem/Davidson/DavidsonSolver.hpp"
#include "ONVBasis/SpinResolvedONVBasis.hpp"
#include "Processing/RDM/RDMCalculator.hpp"
#include "QCMethod/Applications/AtomicDecompositionParameters.hpp"
#include "QCMethod/CI/CI.hpp"
#include "QCMethod/CI/CIEnvironment.hpp"
#include "QCMethod/HF/RHF/DiagonalRHFFockMatrixObjective.hpp"
#include "QCMethod/HF/RHF/RHF.hpp"
#include "QCMethod/HF/RHF/RHFSCFSolver.hpp"
#include "Utilities/units.hpp"


/**
 *  Check if the decomposition of the molecular Hamiltonian for BeH+//STO-3G into atomic contributions works as expected. The dimension of the full spin-resolved Fock space is 441.
 */
BOOST_AUTO_TEST_CASE(decomposition_BeH_cation_STO_3G_Nuclear) {

    // Create the molecular Hamiltonian in an AO basis.
    const GQCP::Nucleus Be {4, 0.0, 0.0, 0.0};
    const GQCP::Nucleus H {1, 0.0, 0.0, GQCP::units::angstrom_to_bohr(1.134)};  // from CCCBDB, STO-3G geometry
    const std::vector<GQCP::Nucleus> nuclei {Be, H};
    const GQCP::Molecule molecule {nuclei, +1};
    const auto N_P = molecule.numberOfElectrons() / 2;

    const GQCP::AtomicDecompositionParameters adp = GQCP::AtomicDecompositionParameters::Nuclear(molecule, "STO-3G");  // the molecular Hamiltonian in its atomic contributions
    const GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};

    auto sq_hamiltonian = adp.molecularHamiltonian();
    const auto K = sq_hamiltonian.dimension();  // number of


    // Transform the molecular Hamiltonian to the canonical RHF basis.
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(molecule.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective {sq_hamiltonian};
    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters();

    const auto& T = rhf_parameters.coefficientMatrix();
    sq_hamiltonian.transform(T);


    // Create the appropriate ONV basis for FCI, specify dense solver and corresponding environment and put them together in the QCMethod to do a FCI calculation.
    const GQCP::SpinResolvedONVBasis onv_basis {K, N_P, N_P};

    const auto initial_guess = GQCP::LinearExpansion<GQCP::SpinResolvedONVBasis>::HartreeFock(onv_basis).coefficients();
    auto environment = GQCP::CIEnvironment::Iterative(sq_hamiltonian, onv_basis, initial_guess);
    auto solver = GQCP::EigenproblemSolver::Davidson();
    const auto qc_structure = GQCP::QCMethod::CI<GQCP::SpinResolvedONVBasis>(onv_basis).optimize(solver, environment);

    const auto electronic_energy = qc_structure.groundStateEnergy();
    const auto& linear_expansion = qc_structure.groundStateParameters();


    // Calculate the RDMs (in the AO basis in which the molecular decomposition parameters are defined) in order to calculate expectation values.
    GQCP::RDMCalculator rdm_calculator {onv_basis};
    rdm_calculator.setCoefficients(linear_expansion.coefficients());

    auto D = rdm_calculator.calculate1RDMs().one_rdm;
    D.basisTransform(T.adjoint());  // T.adjoint() to transform BACK to AO basis

    auto d = rdm_calculator.calculate2RDMs().two_rdm;
    d.basisTransform(T.adjoint());  // T.adjoint() to transform BACK to AO basis


    // Check the decomposed energy values.
    const double repulsion = GQCP::Operator::NuclearRepulsion(molecule).value();

    const double self_energy_a = adp.netAtomic()[0].calculateExpectationValue(D, d);
    const double self_energy_b = adp.netAtomic()[1].calculateExpectationValue(D, d);
    const double interaction_energy_ab = adp.interaction()[0].calculateExpectationValue(D, d) + repulsion;
    const double total_energy_a = adp.atomic()[0].calculateExpectationValue(D, d) + repulsion / 2;
    const double total_energy_b = adp.atomic()[1].calculateExpectationValue(D, d) + repulsion / 2;

    BOOST_CHECK(std::abs(total_energy_a + total_energy_b - electronic_energy - repulsion) < 1.0e-10);
    BOOST_CHECK(std::abs(self_energy_a + self_energy_b + interaction_energy_ab - electronic_energy - repulsion) < 1.0e-10);
    BOOST_CHECK(std::abs(self_energy_a + interaction_energy_ab / 2 - total_energy_a) < 1.0e-10);
    BOOST_CHECK(std::abs(self_energy_b + interaction_energy_ab / 2 - total_energy_b) < 1.0e-10);
}
