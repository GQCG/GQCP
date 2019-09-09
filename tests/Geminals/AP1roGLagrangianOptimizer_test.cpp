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
#define BOOST_TEST_MODULE "AP1roGLagrangianOptimizer"

#include <boost/test/unit_test.hpp>

#include "Geminals/AP1roGLagrangianOptimizer.hpp"
#include "RHF/PlainRHFSCFSolver.hpp"


BOOST_AUTO_TEST_CASE ( h2_631gdp ) {

    // Prepare molecular Hamiltonian in the RHF basis
    auto h2 = GQCP::Molecule::ReadXYZ("data/h2_olsens.xyz");
    GQCP::SingleParticleBasis<double, GQCP::GTOShell> sp_basis (h2, "6-31G**");
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(sp_basis, h2);  // in an AO basis

    GQCP::PlainRHFSCFSolver plain_scf_solver (sq_hamiltonian, h2);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    sq_hamiltonian.transform(rhf.get_C());


    // Optimize the AP1roG Lagrangian, using an initial guess of the geminal coefficients of 0
    GQCP::AP1roGLagrangianOptimizer bivar_solver1 (h2, sq_hamiltonian);
    bivar_solver1.solve();
}
