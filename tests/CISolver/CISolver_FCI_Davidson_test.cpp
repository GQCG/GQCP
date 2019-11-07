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
#define BOOST_TEST_MODULE "DavidsonFCISolver"

#include <boost/test/unit_test.hpp>

#include <iomanip>      // std::setprecision

#include "Basis/transform.hpp"
#include "CISolver/CISolver.hpp"
#include "HamiltonianBuilder/FCI.hpp"
#include "Mathematical/Optimization/DavidsonSolver.hpp"
#include "Operator/SecondQuantized/USQHamiltonian.hpp"
#include "RHF/PlainRHFSCFSolver.hpp"


BOOST_AUTO_TEST_CASE ( FCI_h2_sto3g_dense_vs_Davidson ) {

    // Check if the dense FCI energy is equal to the Davidson (with matvec) FCI energy

    // Create the molecular Hamiltonian in an AO basis
    auto h2 = GQCP::Molecule::ReadXYZ("data/h2_cristina.xyz");
    GQCP::SingleParticleBasis<double, GQCP::GTOShell> sp_basis (h2, "STO-3G");
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(sp_basis, h2);  // in an AO basis
    auto K = sq_hamiltonian.dimension();

    // Create a plain RHF SCF solver and solve the SCF equations
    GQCP::PlainRHFSCFSolver plain_scf_solver (sq_hamiltonian, sp_basis, h2);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    // Transform the Hamiltonian to the RHF orbital basis
    GQCP::basisTransform(sp_basis, sq_hamiltonian, rhf.get_C());

    GQCP::ProductFockSpace fock_space (K, h2.numberOfElectrons()/2, h2.numberOfElectrons()/2);  // dim = 2

    // Create the FCI module
    GQCP::FCI fci (fock_space);
    GQCP::CISolver ci_solver (fci, sq_hamiltonian);

    // Solve Davidson
    GQCP::VectorX<double> initial_g = fock_space.HartreeFockExpansion();
    GQCP::DavidsonSolverOptions davidson_solver_options (initial_g);
    ci_solver.solve(davidson_solver_options);

    // Retrieve the eigenvalues
    auto fci_davidson_eigenvalue = ci_solver.get_eigenpair().get_eigenvalue();

    // Solve Dense
    GQCP::DenseSolverOptions dense_solver_options;
    ci_solver.solve(dense_solver_options);

    // Retrieve the eigenvalues
    auto fci_dense_eigenvalue = ci_solver.get_eigenpair().get_eigenvalue();

    BOOST_CHECK(std::abs(fci_dense_eigenvalue - fci_davidson_eigenvalue) < 1.0e-08);
}


BOOST_AUTO_TEST_CASE ( FCI_H2_6_31Gxx_dense_vs_Davidson ) {

    // Check if the dense FCI energy is equal to the Davidson (with matvec) FCI energy

    // Create the molecular Hamiltonian in an AO basis
    auto h2 = GQCP::Molecule::ReadXYZ("data/h2_cristina.xyz");
    GQCP::SingleParticleBasis<double, GQCP::GTOShell> sp_basis (h2, "6-31G**");
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(sp_basis, h2);  // in an AO basis
    auto K = sq_hamiltonian.dimension();

    // Create a plain RHF SCF solver and solve the SCF equations
    GQCP::PlainRHFSCFSolver plain_scf_solver (sq_hamiltonian, sp_basis, h2);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    // Transform the Hamiltonian to the RHF orbital basis
    GQCP::basisTransform(sp_basis, sq_hamiltonian, rhf.get_C());

    GQCP::ProductFockSpace fock_space (K, h2.numberOfElectrons()/2, h2.numberOfElectrons()/2);  // dim = 100

    // Create the FCI module
    GQCP::FCI fci (fock_space);
    GQCP::CISolver ci_solver (fci, sq_hamiltonian);

    // Solve Davidson
    GQCP::VectorX<double> initial_g = fock_space.HartreeFockExpansion();
    GQCP::DavidsonSolverOptions davidson_solver_options (initial_g);
    ci_solver.solve(davidson_solver_options);

    // Retrieve the eigenvalues
    auto fci_davidson_eigenvalue = ci_solver.get_eigenpair().get_eigenvalue();

    // Solve Dense
    GQCP::DenseSolverOptions dense_solver_options;
    ci_solver.solve(dense_solver_options);

    // Retrieve the eigenvalues
    auto fci_dense_eigenvalue = ci_solver.get_eigenpair().get_eigenvalue();

    BOOST_CHECK(std::abs(fci_dense_eigenvalue - fci_davidson_eigenvalue) < 1.0e-08);
}


BOOST_AUTO_TEST_CASE ( FCI_H2O_STO_3G_dense_vs_Davidson ) {

    // Check if the dense FCI energy is equal to the Davidson (with matvec) FCI energy

    // Create the molecular Hamiltonian in an AO basis
    auto h2o = GQCP::Molecule::ReadXYZ("data/h2o.xyz");
    GQCP::SingleParticleBasis<double, GQCP::GTOShell> sp_basis (h2o, "STO-3G");
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(sp_basis, h2o);  // in an AO basis
    auto K = sq_hamiltonian.dimension();

    // Create a plain RHF SCF solver and solve the SCF equations
    GQCP::PlainRHFSCFSolver plain_scf_solver (sq_hamiltonian, sp_basis, h2o);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    // Transform the Hamiltonian to the RHF orbital basis
    GQCP::basisTransform(sp_basis, sq_hamiltonian, rhf.get_C());

    GQCP::ProductFockSpace fock_space (K, h2o.numberOfElectrons()/2, h2o.numberOfElectrons()/2);  // dim = 441

    // Create the FCI module
    GQCP::FCI fci (fock_space);
    GQCP::CISolver ci_solver (fci, sq_hamiltonian);

    // Solve Davidson
    GQCP::VectorX<double> initial_g = fock_space.HartreeFockExpansion();
    GQCP::DavidsonSolverOptions davidson_solver_options (initial_g);
    ci_solver.solve(davidson_solver_options);

    // Retrieve the eigenvalues
    auto fci_davidson_eigenvalue = ci_solver.get_eigenpair().get_eigenvalue();

    // Solve Dense
    GQCP::DenseSolverOptions dense_solver_options;
    ci_solver.solve(dense_solver_options);

    // Retrieve the eigenvalues
    auto fci_dense_eigenvalue = ci_solver.get_eigenpair().get_eigenvalue();

    BOOST_CHECK(std::abs(fci_dense_eigenvalue - fci_davidson_eigenvalue) < 1.0e-08);
}

BOOST_AUTO_TEST_CASE ( FCI_H6_STO_3G_dense_vs_Davidson ) {

    // Check if the dense FCI energy is equal to the Davidson (with matvec) FCI energy

    // Create the molecular Hamiltonian in an AO basis
    size_t K = 6;
    GQCP::Molecule H6 = GQCP::Molecule::HChain(K, 1.1);
    GQCP::SingleParticleBasis<double, GQCP::GTOShell> sp_basis (H6, "STO-3G");
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(sp_basis, H6);  // in an AO basis

    // Create a plain RHF SCF solver and solve the SCF equations
    GQCP::PlainRHFSCFSolver plain_scf_solver (sq_hamiltonian, sp_basis, H6);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    // Transform the Hamiltonian to the RHF orbital basis
    GQCP::basisTransform(sp_basis, sq_hamiltonian, rhf.get_C());

    GQCP::ProductFockSpace fock_space (K, 3, 2);  // dim = 300

    // Create the FCI module
    GQCP::FCI fci (fock_space);
    GQCP::CISolver ci_solver (fci, sq_hamiltonian);

    // Solve Davidson
    GQCP::VectorX<double> initial_g = fock_space.HartreeFockExpansion();
    GQCP::DavidsonSolverOptions davidson_solver_options (initial_g);
    ci_solver.solve(davidson_solver_options);

    // Retrieve the eigenvalues
    auto fci_davidson_eigenvalue = ci_solver.get_eigenpair().get_eigenvalue();

    // Solve Dense
    GQCP::DenseSolverOptions dense_solver_options;
    ci_solver.solve(dense_solver_options);

    // Retrieve the eigenvalues
    auto fci_dense_eigenvalue = ci_solver.get_eigenpair().get_eigenvalue();

    BOOST_CHECK(std::abs(fci_dense_eigenvalue - fci_davidson_eigenvalue) < 1.0e-08);
}


BOOST_AUTO_TEST_CASE ( FCI_H2O_Unrestricted_Davidson ) {

    // Test if a transformation of a single compoenent for an unrestricted basis results in identical energies for FCI

    // Psi4 and GAMESS' FCI energy (restricted)
    double reference_fci_energy = -75.0129803939602;

    // Create the molecular Hamiltonian in an AO basis
    auto h2o = GQCP::Molecule::ReadXYZ("data/h2o_Psi4_GAMESS.xyz");
    GQCP::SingleParticleBasis<double, GQCP::GTOShell> sp_basis_alpha (h2o, "STO-3G");
    GQCP::SingleParticleBasis<double, GQCP::GTOShell> sp_basis_beta (h2o, "STO-3G");
    auto usq_hamiltonian = GQCP::USQHamiltonian<double>::Molecular(sp_basis_alpha, sp_basis_beta, h2o);  // in an AO basis
    // Create restricted SQHamiltonian to perform RHF
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(sp_basis_alpha, h2o);  // in an AO basis
    auto K = usq_hamiltonian.dimension();

    // Create a plain RHF SCF solver and solve the SCF equations
    GQCP::PlainRHFSCFSolver plain_scf_solver (sq_hamiltonian, sp_basis_alpha, h2o);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    // Transform the Hamiltonian to an orthonormal basis
    GQCP::basisTransform(sp_basis_alpha, sp_basis_beta, usq_hamiltonian, rhf.get_C());

    // Transform the beta component
    // Create stable unitairy matrix
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes (usq_hamiltonian.alphaHamiltonian().core().parameters());
    GQCP::basisTransformBeta(sp_basis_beta, usq_hamiltonian, GQCP::TransformationMatrix<double>(saes.eigenvectors()));


    GQCP::ProductFockSpace fock_space (K, h2o.numberOfElectrons()/2, h2o.numberOfElectrons()/2);  // dim = 441

    // Create the FCI module
    GQCP::FCI fci (fock_space);

    // Davidson solver
    GQCP::DavidsonSolverOptions solver_options(fock_space.HartreeFockExpansion());
    GQCP::VectorX<double> dia = fci.calculateDiagonal(usq_hamiltonian);
    GQCP::VectorFunction matrixVectorProduct = [&fci, &dia, &usq_hamiltonian](const GQCP::VectorX<double>& x) { return fci.matrixVectorProduct(usq_hamiltonian, x, dia); };
    GQCP::DavidsonSolver solver (matrixVectorProduct, dia, solver_options);

    solver.solve();

    // Retrieve the eigenvalues
    auto fci_energy = solver.get_eigenpair().get_eigenvalue();

    // Calculate the total FCI energy
    double internuclear_repulsion_energy = GQCP::Operator::NuclearRepulsion(h2o).value();
    double test_fci_energy = fci_energy + internuclear_repulsion_energy;
    BOOST_CHECK(std::abs(test_fci_energy - (reference_fci_energy)) < 1.0e-06);
}