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
#define BOOST_TEST_MODULE "AP1roGLagrangianNewtonOrbitalOptimizer_test"

#include <boost/test/unit_test.hpp>

#include "Basis/transform.hpp"
#include "RHF/PlainRHFSCFSolver.hpp"
#include "Geminals/AP1roGLagrangianOptimizer.hpp"
#include "Mathematical/Optimization/IterativeIdentitiesHessianModifier.hpp"
#include "OrbitalOptimization/AP1roGLagrangianNewtonOrbitalOptimizer.hpp"


BOOST_AUTO_TEST_CASE ( lih_6_31G_orbital_optimize ) {

    // Construct the molecular Hamiltonian in the RHF basis
    auto lih = GQCP::Molecule::ReadXYZ("data/lih_olsens.xyz");
    GQCP::SingleParticleBasis<double, GQCP::GTOShell> sp_basis (lih, "6-31G");
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(sp_basis, lih);  // in an AO basis

    GQCP::PlainRHFSCFSolver plain_scf_solver (sq_hamiltonian, lih);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();
    basisTransform(sp_basis, sq_hamiltonian, rhf.get_C());

    // Get the initial AP1roG energy
    GQCP::AP1roGLagrangianOptimizer lagrangian_solver (lih, sq_hamiltonian);
    lagrangian_solver.solve();
    double initial_energy = lagrangian_solver.get_electronic_energy();
    auto initial_G = lagrangian_solver.get_geminal_coefficients();

    // Do an AP1roG orbital optimization using Jacobi rotations
    auto hessian_modifier = std::make_shared<GQCP::IterativeIdentitiesHessianModifier>();
    GQCP::AP1roGLagrangianNewtonOrbitalOptimizer orbital_optimizer (initial_G, hessian_modifier, 1.0e-04);
    orbital_optimizer.optimize(sp_basis, sq_hamiltonian);

    double optimized_energy = orbital_optimizer.get_electronic_energy();


    // We don't have reference data, so all we can do is check if orbital optimization lowers the energy
    BOOST_CHECK(optimized_energy < initial_energy);
}
