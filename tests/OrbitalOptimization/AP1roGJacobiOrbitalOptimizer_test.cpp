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
#define BOOST_TEST_MODULE "OO-AP1roG"

#include <boost/test/unit_test.hpp>

#include "Basis/transform.hpp"
#include "Geminals/AP1roG.hpp"
#include "Geminals/AP1roGPSESolver.hpp"
#include "RHF/PlainRHFSCFSolver.hpp"
#include "OrbitalOptimization/AP1roGJacobiOrbitalOptimizer.hpp"


BOOST_AUTO_TEST_CASE ( lih_6_31G_calculateEnergyAfterRotation ) {

    // We have implemented a formula to calculate the rotated AP1roG energy directly, but we have to test it
    // It should be equal to the energy we obtain by rotating the one- and two-electron integrals first


    // Construct the molecular Hamiltonian in the RHF basis
    auto lih = GQCP::Molecule::ReadXYZ("data/lih_olsens.xyz");
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (lih, "6-31G");
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, lih);  // in an AO basis

    GQCP::PlainRHFSCFSolver plain_scf_solver (sq_hamiltonian, spinor_basis, lih);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();
    basisTransform(spinor_basis, sq_hamiltonian, rhf.get_C());


    // Loop over all possible Jacobi pairs for a given (random) angle and check if the analytical result matches the numerical result
    double theta = 56.71;
    size_t K = sq_hamiltonian.dimension();

    for (size_t q = 0; q < K; q++) {  // p and q loop over spatial orbitals
        for (size_t p = q + 1; p < K; p++) {  // p > q
            GQCP::JacobiRotationParameters jacobi_rot_par {p, q, theta};

            GQCP::AP1roGPSESolver pse_solver (lih, sq_hamiltonian);
            pse_solver.solve();
            auto G = pse_solver.get_geminal_coefficients();


            // Calculate the analytical energy after rotation
            GQCP::AP1roGJacobiOrbitalOptimizer orbital_optimizer (G, 1.0e-04);
            orbital_optimizer.calculateJacobiCoefficients(sq_hamiltonian, p, q);
            double E_correction_analytical = orbital_optimizer.calculateScalarFunctionChange(sq_hamiltonian, jacobi_rot_par);


            // Calculate the energy after a numerical rotation (using a Jacobi rotation matrix)
            double E_before = GQCP::calculateAP1roGEnergy(G, sq_hamiltonian);
            auto sq_ham_copy = sq_hamiltonian;
            sq_ham_copy.rotate(jacobi_rot_par);
            double E_after = GQCP::calculateAP1roGEnergy(G, sq_ham_copy);
            double E_correction_numerical = E_after - E_before;


            BOOST_CHECK(std::abs(E_correction_analytical - E_correction_numerical) < 1.0e-08);
        }
    }
}


BOOST_AUTO_TEST_CASE ( lih_6_31G_orbitalOptimize ) {

    // Construct the molecular Hamiltonian in the RHF basis
    auto lih = GQCP::Molecule::ReadXYZ("data/lih_olsens.xyz");
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (lih, "6-31G");
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, lih);  // in an AO basis

    GQCP::PlainRHFSCFSolver plain_scf_solver (sq_hamiltonian, spinor_basis, lih);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();
    basisTransform(spinor_basis, sq_hamiltonian, rhf.get_C());


    // Get the initial AP1roG energy
    GQCP::AP1roGPSESolver pse_solver (lih, sq_hamiltonian);
    pse_solver.solve();
    double initial_energy = pse_solver.get_electronic_energy();
    const auto initial_G = pse_solver.get_geminal_coefficients();


    // Do an AP1roG orbital optimization using Jacobi rotations
    GQCP::AP1roGJacobiOrbitalOptimizer orbital_optimizer (initial_G, 1.0e-04);
    orbital_optimizer.optimize(spinor_basis, sq_hamiltonian);
    double optimized_energy = orbital_optimizer.get_electronic_energy();


    // We don't have reference data, so all we can do is check if orbital optimization lowers the energy
    BOOST_CHECK(optimized_energy < initial_energy);
}
