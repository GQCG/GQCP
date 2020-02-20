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
#define BOOST_TEST_MODULE "properties"

#include <boost/test/unit_test.hpp>

#include "Processing/Properties/properties.hpp"

#include "Basis/Integrals/Interfaces/LibintInterfacer.hpp"
#include "Basis/transform.hpp"
#include "Operator/FirstQuantized/Operator.hpp"
#include "Processing/Properties/RHFElectricalResponseSolver.hpp"
#include "Processing/RDM/RDMCalculator.hpp"
#include "QCMethod/CI/HamiltonianBuilder/FCI.hpp"
#include "QCMethod/HF/DiagonalRHFFockMatrixObjective.hpp"
#include "QCMethod/HF/RHF.hpp"
#include "QCMethod/HF/RHFSCFSolver.hpp"
#include "Utilities/units.hpp"


BOOST_AUTO_TEST_CASE ( dipole_CO_STO_3G ) {

    // Initialize the molecule and molecular Hamiltonian for CO
    GQCP::Nucleus C (6, 0.0, 0.0, 0.0);
    GQCP::Nucleus O (8, 0.0, 0.0, GQCP::units::angstrom_to_bohr(1.145));  // from CCCBDB, STO-3G geometry
    std::vector<GQCP::Nucleus> nuclei {C, O};
    GQCP::Molecule CO (nuclei);

    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (CO, "STO-3G");
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, CO);  // in an AO basis

    size_t K = spinor_basis.numberOfSpatialOrbitals();
    size_t N = CO.numberOfElectrons();

    // Solve the SCF equations
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(CO.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters());
    auto diis_rhf_scf_solver = GQCP::RHFSCFSolver<double>::DIIS();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective (sq_hamiltonian);
    const auto rhf_qc_structure = GQCP::QCMethod::RHF<double>().optimize(objective, diis_rhf_scf_solver, rhf_environment);
    const auto rhf_parameters = rhf_qc_structure.groundStateParameters();

    double total_energy = rhf_qc_structure.groundStateEnergy() + GQCP::Operator::NuclearRepulsion(CO).value();
    BOOST_REQUIRE(std::abs(total_energy - (-111.225)) < 1.0e-02);  // from CCCBDB, require a correct RHF solution to be found


    // Calculate the RHF 1-RDM in MO basis
    auto D = GQCP::QCModel::RHF<double>::calculateOrthonormalBasis1RDM(K, N);
    auto D_AO = GQCP::QCModel::RHF<double>::calculateScalarBasis1RDM(rhf_parameters.coefficientMatrix(), N);

    // Calculate the dipole integrals, and transform them to the MO basis
    auto dipole_op = spinor_basis.quantize(GQCP::Operator::ElectronicDipole());
    dipole_op.transform(rhf_parameters.coefficientMatrix());

    GQCP::Vector<double, 3> total_dipole_moment = GQCP::Operator::NuclearDipole(CO).value() + dipole_op.calculateExpectationValue(D);
    BOOST_CHECK(std::abs(total_dipole_moment.norm() - (0.049)) < 1.0e-03);
}


BOOST_AUTO_TEST_CASE ( dipole_N2_STO_3G ) {

    // Check that the dipole moment of N2 is zero

    // Initialize the molecule and the molecular Hamiltonian for N2
    GQCP::Nucleus N_1 (7, 0.0, 0.0, 0.0);
    GQCP::Nucleus N_2 (7, 0.0, 0.0, GQCP::units::angstrom_to_bohr(1.134));  // from CCCBDB, STO-3G geometry
    GQCP::Molecule N2 ({N_1, N_2});

    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (N2, "STO-3G");
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, N2);  // in an AO basis

    size_t K = spinor_basis.numberOfSpatialOrbitals();
    size_t N = N2.numberOfElectrons();

    // Solve the SCF equations
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(N2.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective (sq_hamiltonian);
    const auto rhf_qc_structure = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment);
    const auto rhf_parameters = rhf_qc_structure.groundStateParameters();

    double total_energy = rhf_qc_structure.groundStateEnergy() + GQCP::Operator::NuclearRepulsion(N2).value();
    BOOST_REQUIRE(std::abs(total_energy - (-107.500654)) < 1.0e-05);  // from CCCBDB, require a correct RHF solution to be found


    // Calculate the RHF 1-RDM in MO basis
    auto D = GQCP::QCModel::RHF<double>::calculateOrthonormalBasis1RDM(K, N);
    auto D_AO = GQCP::QCModel::RHF<double>::calculateScalarBasis1RDM(rhf_parameters.coefficientMatrix(), N);

    // Calculate the dipole integrals, and transform them to the MO basis
    auto dipole_op = spinor_basis.quantize(GQCP::Operator::ElectronicDipole());
    dipole_op.transform(rhf_parameters.coefficientMatrix());

    GQCP::Vector<double, 3> total_dipole_moment = GQCP::Operator::NuclearDipole(N2).value() + dipole_op.calculateExpectationValue(D);
    BOOST_CHECK(std::abs(total_dipole_moment.norm() - (0.0)) < 1.0e-08);
}


/**
 *  Check the calculation of the zz-component of the polarizability for H2 with a reference value from Psi4-numpy
 * 
 *  Note that the reference value is generated from Psi4-numpy, with a fix for the Fockian matrix
 */
BOOST_AUTO_TEST_CASE ( h2_polarizability_RHF ) {

    // Initialize the reference value
    const double ref_alpha_zz = 1.08428;


    // Initialize the molecule and the Hamiltonian in the AO basis
    GQCP::Nucleus H1 (1, 0.0, 0.0,  0.0);
    GQCP::Nucleus H2 (1, 0.0, 0.0,  0.5);
    GQCP::Molecule h2 ({H1, H2}, 0);

    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (h2, "STO-3G");
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, h2);  // in the AO basis


    // Do the RHF calculation to get the canonical RHF orbitals
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(h2.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective (sq_hamiltonian);
    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters();


    // Transform the orbitals to the RHF basis and prepare the dipole integrals in the RHF basis
    GQCP::basisTransform(spinor_basis, sq_hamiltonian, rhf_parameters.coefficientMatrix());
    const auto dipole_op = spinor_basis.quantize(GQCP::Operator::ElectronicDipole());


    // Find the RHF wave function response
    GQCP::RHFElectricalResponseSolver cphf_solver (h2.numberOfElectrons()/2);
    const auto x = cphf_solver.calculateWaveFunctionResponse(sq_hamiltonian, dipole_op);


    // Calculate the RHF polarizability
    const auto F_p = cphf_solver.calculateParameterResponseForce(dipole_op);
    const auto alpha = GQCP::calculateElectricPolarizability(F_p, x);
    const auto alpha_zz = alpha(2,2);


    BOOST_CHECK(std::abs(alpha_zz - ref_alpha_zz) < 1.0e-05);
}


/**
 *  Test the Dyson algorithm against manually calculated coefficients for two (normalized) toy wave functions
 */
BOOST_AUTO_TEST_CASE ( dyson_coefficients ) {

    // Set up the manually calculated references
    GQCP::VectorX<double> reference_amplitudes_beta = GQCP::VectorX<double>::Zero(2); 
    reference_amplitudes_beta << 0.537653264399, 0.794791398869;

    GQCP::VectorX<double> reference_amplitudes_alpha = GQCP::VectorX<double>::Zero(2); 
    reference_amplitudes_alpha << 0.39739531532399996, 0.9116729926689999;

    // Set up the toy wave functions
    const size_t K = 2;
    const size_t N = 2;
 
    const GQCP::SpinResolvedONVBasis fock_space1 (K, N/2, N/2);
    const GQCP::SpinResolvedONVBasis fock_space2 (K, N/2, N/2-1);
    const GQCP::SpinResolvedONVBasis fock_space3 (K, N/2-1, N/2);
  
    GQCP::VectorX<double> coeffs1 = GQCP::VectorX<double>::Zero(4); 
    coeffs1 << 0.182574, 0.365148, 0.547723, 0.730297;
    GQCP::VectorX<double> coeffs2 = GQCP::VectorX<double>::Zero(2); 
    coeffs2 << 0.640184, 0.768221;
    
    const auto linear_expansion1 = GQCP::LinearExpansion<GQCP::SpinResolvedONVBasis>(fock_space1, coeffs1);
    const auto linear_expansion2 = GQCP::LinearExpansion<GQCP::SpinResolvedONVBasis>(fock_space2, coeffs2);
    const auto linear_expansion3 = GQCP::LinearExpansion<GQCP::SpinResolvedONVBasis>(fock_space3, coeffs2);

    // Calculate the coefficients of the Dyson orbitals and check with the reference
    const auto dyson_coefficients_beta = GQCP::calculateDysonOrbitalCoefficients(linear_expansion1, linear_expansion2);  // coefficients with a difference in beta occupation
    const auto dyson_coefficients_alpha = GQCP::calculateDysonOrbitalCoefficients(linear_expansion1, linear_expansion3);  // coefficients with a difference in alpha occupation

    BOOST_CHECK(dyson_coefficients_beta.isApprox(reference_amplitudes_beta, 1.0e-6));
    BOOST_CHECK(dyson_coefficients_alpha.isApprox(reference_amplitudes_alpha, 1.0e-6));
}
