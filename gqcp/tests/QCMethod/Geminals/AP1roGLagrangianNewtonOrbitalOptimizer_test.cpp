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

#define BOOST_TEST_MODULE "AP1roGLagrangianNewtonOrbitalOptimizer_test"

#include <boost/test/unit_test.hpp>

#include "Basis/Transformations/transform.hpp"
#include "Mathematical/Optimization/Minimization/IterativeIdentitiesHessianModifier.hpp"
#include "Mathematical/Optimization/NonLinearEquation/NonLinearEquationSolver.hpp"
#include "QCMethod/Geminals/AP1roG.hpp"
#include "QCMethod/Geminals/AP1roGLagrangianNewtonOrbitalOptimizer.hpp"
#include "QCMethod/Geminals/PSEnvironment.hpp"
#include "QCMethod/HF/RHF/DiagonalRHFFockMatrixObjective.hpp"
#include "QCMethod/HF/RHF/RHF.hpp"
#include "QCMethod/HF/RHF/RHFSCFSolver.hpp"
#include "QCModel/Geminals/AP1roG.hpp"


/**
 *  Since we don't have OO reference data, all we can do is check if the orbital optimization lowers the energy.
 */
BOOST_AUTO_TEST_CASE(lih_6_31G_orbital_optimize) {

    // Construct the molecular Hamiltonian in the RHF basis
    const auto lih = GQCP::Molecule::ReadXYZ("data/lih_olsens.xyz");
    const auto N_P = lih.numberOfElectrons() / 2;
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis {lih, "6-31G"};
    auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spinor_basis, lih);  // in an AO basis

    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(lih.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective {sq_hamiltonian};
    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters();
    basisTransform(spinor_basis, sq_hamiltonian, rhf_parameters.coefficientMatrix());


    // Get the initial AP1roG solution.
    auto solver = GQCP::NonLinearEquationSolver<double>::Newton();
    auto environment = GQCP::PSEnvironment::AP1roG(sq_hamiltonian, N_P);  // zero initial guess
    const auto qc_structure = GQCP::QCMethod::AP1roG(sq_hamiltonian, N_P).optimize(solver, environment);

    const auto G_initial = qc_structure.groundStateParameters().geminalCoefficients();
    const auto initial_energy = qc_structure.groundStateEnergy();


    // Do an AP1roG orbital optimization using a Newton-based algorithm.
    auto hessian_modifier = std::make_shared<GQCP::IterativeIdentitiesHessianModifier>();
    GQCP::AP1roGLagrangianNewtonOrbitalOptimizer orbital_optimizer {G_initial, hessian_modifier, 1.0e-04};
    orbital_optimizer.optimize(spinor_basis, sq_hamiltonian);
    const auto optimized_energy = orbital_optimizer.electronicEnergy();

    // We don't have reference data, so all we can do is check if orbital optimization lowers the energy
    BOOST_CHECK(optimized_energy < initial_energy);
}
