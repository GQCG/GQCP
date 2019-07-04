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
#define BOOST_TEST_MODULE "WaveFunction"

#include <boost/test/unit_test.hpp>

#include "WaveFunction/WaveFunction.hpp"
#include "FockSpace/FockSpace.hpp"
#include "CISolver/CISolver.hpp"
#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "HamiltonianBuilder/FCI.hpp"
#include "Molecule.hpp"
#include "RHF/PlainRHFSCFSolver.hpp"

BOOST_AUTO_TEST_CASE ( shannon_entropy ) {

    // Set up a test Fock space
    GQCP::FockSpace fock_space (8, 3);  // K = 8, N = 3


    // Check the entropy of a Hartree-Fock expansion
    GQCP::WaveFunction hartree_fock_expansion (fock_space, fock_space.HartreeFockExpansion());
    BOOST_CHECK(hartree_fock_expansion.calculateShannonEntropy() < 1.0e-12);  // should be 0


    // Check the maximal entropy (corresponding to a wave function with all equal coefficients different from zero)
    GQCP::WaveFunction constant_expansion (fock_space, fock_space.constantExpansion());
    double reference_entropy = std::log2(fock_space.get_dimension());  // manual derivation
    BOOST_CHECK(std::abs(constant_expansion.calculateShannonEntropy() - reference_entropy) < 1.0e-12);
}

BOOST_AUTO_TEST_CASE ( transform_wave_function ) {

    // Produce a wave function, transform it then pair it against second produced wave function from a transformed basis.
    std::cout<<"Checkpoint0";

        // Create a molecule
    GQCP::Molecule hchain = GQCP::Molecule::HChain(3, 0.742, -1);

    // Create the molecular Hamiltonian parameters for this molecule and basis
    auto mol_ham_par = GQCP::HamiltonianParameters<double>::Molecular(hchain, "STO-3G");
    auto K = mol_ham_par.get_K();
    auto N_P = hchain.get_N()/2;

    std::cout<<"K :"<<K<<" ";
    std::cout<<"N_P :"<<N_P<<" ";

    // Create a plain RHF SCF solver and solve the SCF equations
    GQCP::PlainRHFSCFSolver plain_scf_solver (mol_ham_par, hchain);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    std::cout<<"Checkpoint";

    mol_ham_par.LowdinOrthonormalize();
    GQCP::ProductFockSpace fock_space (K, N_P, N_P);
    GQCP::FCI fci (fock_space);

    GQCP::DenseSolverOptions solver_options;
    GQCP::CISolver ci_solver (fci, mol_ham_par);
    ci_solver.solve(solver_options);

    // Retrieve the wave function and transform it
    auto wavefunction1 = ci_solver.makeWavefunction();
    std::cout<<std::endl<<wavefunction1.get_coefficients();
    std::cout<<std::endl<<"-------------"<<std::endl;
    GQCP::SquareMatrix<double> U_random = GQCP::SquareMatrix<double>::RandomUnitary(K);
    std::cout<<U_random;
    std::cout<<"Checkpoint2";


    wavefunction1.basisTransform(U_random);

    GQCP::SquareMatrix<double> U_T = GQCP::SquareMatrix<double>(U_random.transpose());
    // Generate a new wave function by rotating the basis and performing the FCI again.
    mol_ham_par.rotate(U_random);
    GQCP::CISolver ci_solver2 (fci, mol_ham_par);
    ci_solver2.solve(solver_options);

    auto wavefunction2 = ci_solver2.makeWavefunction();

    std::cout<<std::endl<<wavefunction2.get_coefficients();
    std::cout<<std::endl<<"-------------"<<std::endl;
    std::cout<<std::endl<<wavefunction1.get_coefficients();

    // Check if they deviate
    BOOST_CHECK(wavefunction2.get_coefficients().isApprox(wavefunction1.get_coefficients()));



    BOOST_CHECK(true);
}
