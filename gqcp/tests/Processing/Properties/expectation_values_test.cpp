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
#define BOOST_TEST_MODULE "expectation_values"

#include <boost/test/unit_test.hpp>

#include "Processing/Properties/expectation_values.hpp"

#include "Basis/Integrals/Interfaces/LibintInterfacer.hpp"
#include "Basis/transform.hpp"
#include "Mathematical/Optimization/Eigenproblem/DavidsonSolver.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "Operator/SecondQuantized/USQHamiltonian.hpp"
#include "Processing/RDM/RDMCalculator.hpp"
#include "QCMethod/CI/CISolver.hpp"
#include "QCMethod/CI/HamiltonianBuilder/DOCI.hpp"
#include "QCMethod/CI/HamiltonianBuilder/FCI.hpp"
#include "QCMethod/HF/DiagonalRHFFockMatrixObjective.hpp"
#include "QCMethod/HF/RHF.hpp"
#include "QCMethod/HF/RHFSCFSolver.hpp"
#include "Utilities/units.hpp"


BOOST_AUTO_TEST_CASE ( mulliken_N2_STO_3G ) {

    // Check that the mulliken population of N2 is 14 (N)

    // Initialize the molecule and the molecular Hamiltonian for N2
    GQCP::Nucleus N_1 (7, 0.0, 0.0, 0.0);
    GQCP::Nucleus N_2 (7, 0.0, 0.0, GQCP::units::angstrom_to_bohr(1.134));  // from CCCBDB, STO-3G geometry
    std::vector<GQCP::Nucleus> nuclei {N_1, N_2};
    GQCP::Molecule N2 (nuclei);

    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (N2, "STO-3G");
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, N2);  // in an AO basis
    size_t K = sq_hamiltonian.dimension();

    // We include all basis functions
    std::vector<size_t> gto_list (K);
    for(size_t i = 0; i<K; i++){
        gto_list[i] = i;
    }

    GQCP::ScalarSQOneElectronOperator<double> mulliken = spinor_basis.calculateMullikenOperator(gto_list);

    size_t N = N2.numberOfElectrons();

    // Create a 1-RDM for N2
    GQCP::OneRDM<double> one_rdm = GQCP::QCModel::RHF<double>::calculateOrthonormalBasis1RDM(K, N);

    double mulliken_population = mulliken.calculateExpectationValue(one_rdm)(0);
    BOOST_CHECK(std::abs(mulliken_population - (N)) < 1.0e-06);


    // Repeat this for a DOCI-RDM

    // Solve the SCF equations
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(N2.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective (sq_hamiltonian);
    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters();

    sq_hamiltonian.transform(rhf_parameters.coefficientMatrix());

    GQCP::SpinUnresolvedONVBasis onv_basis (K, N/2);
    GQCP::DOCI doci (onv_basis);

    GQCP::CISolver ci_solver (doci, sq_hamiltonian);

    GQCP::DenseSolverOptions solver_options;
    ci_solver.solve(solver_options);

    GQCP::RDMCalculator rdm_calculator (ci_solver.makeLinearExpansion());

    GQCP::OneRDMs<double> one_rdms = rdm_calculator.calculate1RDMs();

    double mulliken_population_2 = mulliken.calculateExpectationValue(one_rdms.one_rdm)(0);
    BOOST_CHECK(std::abs(mulliken_population_2 - (N)) < 1.0e-06);
}


/*
 *  Perturb the S_z for a single atom in the diatomic molecule NO+
 *  After the perturbation, the total S_z for the diatomic molecule should still be 0 while the fragment S_z should have changed
 */ 
BOOST_AUTO_TEST_CASE ( S_z_constrained_NOplus_STO_3G ) {

    // Initialize the molecule and the molecular Hamiltonian for NO+
    GQCP::Nucleus N (7, 0.0, 0.0, 0.0);
    GQCP::Nucleus O (8, 0.0, 0.0, 2);  
    std::vector<GQCP::Nucleus> nuclei {N, O};
    GQCP::Molecule NOplus (nuclei, +1);

    GQCP::USpinorBasis<double, GQCP::GTOShell> uspinor_basis (NOplus, "STO-3G");
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (NOplus, "STO-3G");
    auto usq_hamiltonian = GQCP::USQHamiltonian<double>::Molecular(uspinor_basis, NOplus);  // in an AO basis
    
    // Create restricted Hamiltonian to perform RHF
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, NOplus);  // in an AO basis
    size_t K = sq_hamiltonian.dimension();

    // Basis functions for O
    std::vector<size_t> gto_O (K/2);
    for(size_t i = K/2; i<K; i++){
        gto_O[i-K/2] = i;
    }

    // Basis functions for N
    std::vector<size_t> gto_N (K/2);
    for(size_t i = 0; i<K/2; i++){
        gto_N[i] = i;
    }
    
    size_t Ne = NOplus.numberOfElectrons();
    
    // Solve the SCF equations
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(NOplus.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters());
    auto diis_rhf_scf_solver = GQCP::RHFSCFSolver<double>::DIIS();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective (sq_hamiltonian);
    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, diis_rhf_scf_solver, rhf_environment).groundStateParameters();

    GQCP::basisTransform(uspinor_basis, usq_hamiltonian, rhf_parameters.coefficientMatrix());

    // Calculate the atomic spin-z operator
    GQCP::ScalarSQOneElectronOperator<double> sq_N_Sz_alpha = uspinor_basis.calculateAtomicSpinZ(gto_N, GQCP::SpinComponent::ALPHA);
    GQCP::ScalarSQOneElectronOperator<double> sq_N_Sz_beta = uspinor_basis.calculateAtomicSpinZ(gto_N, GQCP::SpinComponent::BETA);
    GQCP::ScalarSQOneElectronOperator<double> sq_O_Sz_alpha = uspinor_basis.calculateAtomicSpinZ(gto_O, GQCP::SpinComponent::ALPHA);
    GQCP::ScalarSQOneElectronOperator<double> sq_O_Sz_beta = uspinor_basis.calculateAtomicSpinZ(gto_O, GQCP::SpinComponent::BETA);

    // Create a constraint for the spin-z on N
    auto constrained = usq_hamiltonian.constrain(sq_N_Sz_alpha, 0.5, GQCP::SpinComponent::ALPHA);
    constrained = constrained.constrain(sq_N_Sz_beta, 0.5, GQCP::SpinComponent::BETA);

    // Do the FCI calculation 
    GQCP::SpinResolvedONVBasis onv_basis (K, Ne/2, Ne/2);
    
    GQCP::DavidsonSolverOptions solver_options (onv_basis.hartreeFockExpansion());
    GQCP::VectorX<double> dia = onv_basis.evaluateOperatorDiagonal(constrained);
    GQCP::VectorFunction<double> matrixVectorProduct = [&onv_basis, &dia, &constrained](const GQCP::VectorX<double>& x) { return onv_basis.evaluateOperatorMatrixVectorProduct(constrained, x, dia); };
    GQCP::DavidsonSolver ds_solver (matrixVectorProduct, dia, solver_options);

    ds_solver.solve();

    // Calculate the RDMs in order to evaluate expectation values
    GQCP::RDMCalculator rdm_calc(onv_basis);
    rdm_calc.set_coefficients(ds_solver.get_eigenpair().get_eigenvector());

    auto one_rdms = rdm_calc.calculate1RDMs();

    // Calculate the spin density matrix
    GQCP::OneRDM<double> spin_d = one_rdms.spinDensityRDM();

    // Evaluate S_z for O and N
    double N_Sz = sq_N_Sz_alpha.calculateExpectationValue(spin_d)[0];
    double O_Sz = sq_O_Sz_alpha.calculateExpectationValue(spin_d)[0];
    
    // Check that the total Sz (N_Sz + O_Sz) still equals 0
    BOOST_CHECK(std::abs(N_Sz + O_Sz) < 1.0e-06);

    // Check that the Sz for a single fragment has changed
    BOOST_CHECK(std::abs(N_Sz) > 1.0e-06);
}


/*
 *  Calculate Sz and S^2 values for O2 in a restricted SpinUnresolvedONV basis (alpha == beta)
 */ 
BOOST_AUTO_TEST_CASE ( spin_O2 ) {

    // Initialize the molecule and the molecular Hamiltonian for O2
    GQCP::Nucleus O_1 (8, 0.0, 0.0, 0.0);
    GQCP::Nucleus O_2 (8, 0.0, 0.0, 2);
    std::vector<GQCP::Nucleus> nuclei {O_1, O_2};
    GQCP::Molecule O2 (nuclei);

    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (O2, "STO-3G");
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, O2);  // in an AO basis
    size_t K = sq_hamiltonian.dimension();
    size_t N = O2.numberOfElectrons();

    // Solve the SCF equations
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(O2.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters());
    auto diis_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();  // the DIIS SCF solver seems to find a wrong minimum, so use a plain solver instead
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective (sq_hamiltonian);
    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, diis_rhf_scf_solver, rhf_environment).groundStateParameters();

    sq_hamiltonian.transform(rhf_parameters.coefficientMatrix());

    GQCP::SpinResolvedONVBasis onv_basis (K, N/2, N/2);
    GQCP::FCI fci (onv_basis);
 
    GQCP::CISolver ci_solver (fci, sq_hamiltonian);

    GQCP::DenseSolverOptions solver_options; // Dense is required, Davidson will not converge to the lowest eigenstate
    ci_solver.solve(solver_options);
    
    GQCP::RDMCalculator rdm_calculator (ci_solver.makeLinearExpansion());
    
    GQCP::OneRDMs<double> one_rdms = rdm_calculator.calculate1RDMs();
    GQCP::TwoRDMs<double> two_rdms = rdm_calculator.calculate2RDMs();
    
    double s_squared = GQCP::calculateSpinSquared<double>(one_rdms, two_rdms);
    double s_z = GQCP::calculateSpinZ<double>(one_rdms);

    // <S^2> should be 2 (S=1) because the ground state for O2 is a biradical triplet.
    BOOST_CHECK(std::abs(s_squared - 2) < 1.0e-06);

    // In the restricted SpinUnresolvedONV basis, alpha = beta, hence the expectation value of the z-component of the spin operator should be zero
    BOOST_CHECK(std::abs(s_z - 0) < 1.0e-06);
}
