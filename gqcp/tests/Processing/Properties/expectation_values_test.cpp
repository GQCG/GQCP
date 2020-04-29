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

#define BOOST_TEST_MODULE "expectation_values"

#include <boost/test/unit_test.hpp>

#include "Basis/Integrals/Interfaces/LibintInterfacer.hpp"
#include "Basis/transform.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "Operator/SecondQuantized/USQHamiltonian.hpp"
#include "Processing/Properties/expectation_values.hpp"
#include "Processing/RDM/RDMCalculator.hpp"
#include "QCMethod/CI/HamiltonianBuilder/DOCI.hpp"
#include "QCMethod/CI/HamiltonianBuilder/FCI.hpp"
#include "QCMethod/HF/DiagonalRHFFockMatrixObjective.hpp"
#include "QCMethod/HF/RHF.hpp"
#include "QCMethod/HF/RHFSCFSolver.hpp"
#include "Utilities/units.hpp"


/**
 *  Check that the total RHF Mulliken population of N2 is 14.
 */
BOOST_AUTO_TEST_CASE(mulliken_N2_STO_3G) {

    // Initialize the molecular Hamiltonian for N2 in the Löwdin-orthonormalized basis.
    const GQCP::Nucleus N1 {7, 0.0, 0.0, 0.0};
    const GQCP::Nucleus N2 {7, 0.0, 0.0, GQCP::units::angstrom_to_bohr(1.134)};  // from CCCBDB, STO-3G geometry
    const GQCP::Molecule molecule {{N1, N2}};
    const auto N = molecule.numberOfElectrons();

    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    spinor_basis.lowdinOrthonormalize();
    const auto K = spinor_basis.numberOfSpatialOrbitals();

    const auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, molecule);


    // Calculate the Mulliken population operator in this spinor basis.

    // To calculate the total Mulliken population operator, we have to include all basis functions.
    std::vector<size_t> basis_functions;
    basis_functions.reserve(K);
    for (size_t i = 0; i < K; i++) {
        basis_functions.push_back(i);
    }
    const auto mulliken_op = spinor_basis.calculateMullikenOperator(basis_functions);  // 'op' for 'operator'


    // Create the RHF 1-DM for N2 and check the total Mulliken operator.
    const auto D = GQCP::QCModel::RHF<double>::calculateOrthonormalBasis1RDM(K, N);

    const auto mulliken_population = mulliken_op.calculateExpectationValue(D)(0);
    BOOST_CHECK(std::abs(mulliken_population - N) < 1.0e-12);
}


// /*
//  *  Perturb the S_z for a single atom in the diatomic molecule NO+
//  *  After the perturbation, the total S_z for the diatomic molecule should still be 0 while the fragment S_z should have changed
//  */
// BOOST_AUTO_TEST_CASE ( S_z_constrained_NOplus_STO_3G ) {

//     // Initialize the molecule and the molecular Hamiltonian for NO+
//     GQCP::Nucleus N (7, 0.0, 0.0, 0.0);
//     GQCP::Nucleus O (8, 0.0, 0.0, 2);
//     std::vector<GQCP::Nucleus> nuclei {N, O};
//     GQCP::Molecule NOplus (nuclei, +1);

//     GQCP::USpinorBasis<double, GQCP::GTOShell> uspinor_basis (NOplus, "STO-3G");
//     GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (NOplus, "STO-3G");
//     auto usq_hamiltonian = GQCP::USQHamiltonian<double>::Molecular(uspinor_basis, NOplus);  // in an AO basis

//     // Create restricted Hamiltonian to perform RHF
//     auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, NOplus);  // in an AO basis
//     size_t K = sq_hamiltonian.dimension();

//     // Basis functions for O
//     std::vector<size_t> gto_O (K/2);
//     for(size_t i = K/2; i<K; i++){
//         gto_O[i-K/2] = i;
//     }

//     // Basis functions for N
//     std::vector<size_t> gto_N (K/2);
//     for(size_t i = 0; i<K/2; i++){
//         gto_N[i] = i;
//     }

//     size_t Ne = NOplus.numberOfElectrons();

//     // Solve the SCF equations
//     auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(NOplus.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters());
//     auto diis_rhf_scf_solver = GQCP::RHFSCFSolver<double>::DIIS();
//     const GQCP::DiagonalRHFFockMatrixObjective<double> objective (sq_hamiltonian);
//     const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, diis_rhf_scf_solver, rhf_environment).groundStateParameters();

//     GQCP::basisTransform(uspinor_basis, usq_hamiltonian, rhf_parameters.coefficientMatrix());

//     // Calculate the atomic spin-z operator
//     GQCP::ScalarSQOneElectronOperator<double> sq_N_Sz_alpha = uspinor_basis.calculateAtomicSpinZ(gto_N, GQCP::SpinComponent::ALPHA);
//     GQCP::ScalarSQOneElectronOperator<double> sq_N_Sz_beta = uspinor_basis.calculateAtomicSpinZ(gto_N, GQCP::SpinComponent::BETA);
//     GQCP::ScalarSQOneElectronOperator<double> sq_O_Sz_alpha = uspinor_basis.calculateAtomicSpinZ(gto_O, GQCP::SpinComponent::ALPHA);
//     GQCP::ScalarSQOneElectronOperator<double> sq_O_Sz_beta = uspinor_basis.calculateAtomicSpinZ(gto_O, GQCP::SpinComponent::BETA);

//     // Create a constraint for the spin-z on N
//     auto constrained = usq_hamiltonian.constrain(sq_N_Sz_alpha, 0.5, GQCP::SpinComponent::ALPHA);
//     constrained = constrained.constrain(sq_N_Sz_beta, 0.5, GQCP::SpinComponent::BETA);

//     // Do the FCI calculation
//     GQCP::SpinResolvedONVBasis onv_basis (K, Ne/2, Ne/2);

//     GQCP::DavidsonSolverOptions solver_options (onv_basis.hartreeFockExpansion());
//     GQCP::VectorX<double> dia = onv_basis.evaluateOperatorDiagonal(constrained);
//     GQCP::VectorFunction<double> matrixVectorProduct = [&onv_basis, &dia, &constrained](const GQCP::VectorX<double>& x) { return onv_basis.evaluateOperatorMatrixVectorProduct(constrained, x, dia); };
//     GQCP::DavidsonSolver ds_solver (matrixVectorProduct, dia, solver_options);

//     ds_solver.solve();

//     // Calculate the RDMs in order to evaluate expectation values
//     GQCP::RDMCalculator rdm_calc(onv_basis);
//     rdm_calc.set_coefficients(ds_solver.get_eigenpair().get_eigenvector());

//     auto one_rdms = rdm_calc.calculate1RDMs();

//     // Calculate the spin density matrix
//     GQCP::OneRDM<double> spin_d = one_rdms.spinDensityRDM();

//     // Evaluate S_z for O and N
//     double N_Sz = sq_N_Sz_alpha.calculateExpectationValue(spin_d)[0];
//     double O_Sz = sq_O_Sz_alpha.calculateExpectationValue(spin_d)[0];

//     // Check that the total Sz (N_Sz + O_Sz) still equals 0
//     BOOST_CHECK(std::abs(N_Sz + O_Sz) < 1.0e-06);

//     // Check that the Sz for a single fragment has changed
//     BOOST_CHECK(std::abs(N_Sz) > 1.0e-06);
// }


/*
 *  Calculate Sz and S^2 values for O2 in a restricted SpinUnresolvedONV basis (alpha == beta)
 */
// BOOST_AUTO_TEST_CASE ( spin_operators_O2 ) {

//     // Initialize the molecule and the molecular Hamiltonian for O2
//     GQCP::Nucleus O_1 (8, 0.0, 0.0, 0.0);
//     GQCP::Nucleus O_2 (8, 0.0, 0.0, 2);
//     GQCP::Molecule molecule (nuclei);

//     GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (O2, "STO-3G");
//     auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, O2);  // in an AO basis
//     size_t K = sq_hamiltonian.dimension();
//     size_t N = O2.numberOfElectrons();

//     // Solve the SCF equations
//     auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(O2.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters());
//     auto diis_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();  // the DIIS SCF solver seems to find a wrong minimum, so use a plain solver instead
//     const GQCP::DiagonalRHFFockMatrixObjective<double> objective (sq_hamiltonian);
//     const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, diis_rhf_scf_solver, rhf_environment).groundStateParameters();

//     sq_hamiltonian.transform(rhf_parameters.coefficientMatrix());

//     GQCP::SpinResolvedONVBasis onv_basis (K, N/2, N/2);
//     GQCP::FCI fci (onv_basis);

//     GQCP::CISolver ci_solver (fci, sq_hamiltonian);

//     GQCP::DenseSolverOptions solver_options; // Dense is required, Davidson will not converge to the lowest eigenstate
//     ci_solver.solve(solver_options);

//     GQCP::RDMCalculator rdm_calculator (ci_solver.makeLinearExpansion());

//     GQCP::OneRDMs<double> one_rdms = rdm_calculator.calculate1RDMs();
//     GQCP::TwoRDMs<double> two_rdms = rdm_calculator.calculate2RDMs();

//     double s_squared = GQCP::calculateSpinSquared<double>(one_rdms, two_rdms);
//     double s_z = GQCP::calculateSpinZ<double>(one_rdms);

//     // <S^2> should be 2 (S=1) because the ground state for O2 is a biradical triplet.
//     BOOST_CHECK(std::abs(s_squared - 2) < 1.0e-06);

//     // In the restricted SpinUnresolvedONV basis, alpha = beta, hence the expectation value of the z-component of the spin operator should be zero
//     BOOST_CHECK(std::abs(s_z - 0) < 1.0e-06);
// }
