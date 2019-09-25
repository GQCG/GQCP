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

#include "Basis/Integrals/Interfaces/LibintInterfacer.hpp"
#include "CISolver/CISolver.hpp"
#include "HamiltonianBuilder/FCI.hpp"
#include "Operator/FirstQuantized/Operator.hpp"
#include "Properties/properties.hpp"
#include "RHF/DIISRHFSCFSolver.hpp"
#include "RHF/PlainRHFSCFSolver.hpp"
#include "RDM/RDMCalculator.hpp"
#include "units.hpp"


BOOST_AUTO_TEST_CASE ( dipole_CO_STO_3G ) {

    // Initialize the molecule and molecular Hamiltonian parameters for CO
    GQCP::Nucleus C (6, 0.0, 0.0, 0.0);
    GQCP::Nucleus O (8, 0.0, 0.0, GQCP::units::angstrom_to_bohr(1.145));  // from CCCBDB, STO-3G geometry
    std::vector<GQCP::Nucleus> nuclei {C, O};
    GQCP::Molecule CO (nuclei);

    auto ao_basis = std::make_shared<GQCP::ScalarBasis<GQCP::GTOShell>>(CO, "STO-3G");
    auto ao_mol_ham_par = GQCP::HamiltonianParameters<double>::Molecular(ao_basis);

    size_t K = ao_basis->numberOfBasisFunctions();
    size_t N = CO.numberOfElectrons();

    // Solve the SCF equations
    GQCP::DIISRHFSCFSolver diis_scf_solver (ao_mol_ham_par, CO);
    diis_scf_solver.solve();
    auto rhf = diis_scf_solver.get_solution();

    double total_energy = rhf.get_electronic_energy() + GQCP::Operator::NuclearRepulsion(CO).value();
    BOOST_REQUIRE(std::abs(total_energy - (-111.225)) < 1.0e-02);  // from CCCBDB, require a correct RHF solution to be found


    // Calculate the RHF 1-RDM in MO basis
    auto D = GQCP::calculateRHF1RDM(K, N);
    auto D_AO = GQCP::calculateRHFAO1RDM(rhf.get_C(), N);

    // Calculate the dipole integrals, and transform them to the MO basis
    const GQCP::SingleParticleBasis<double, GQCP::GTOShell> sp_basis (*ao_basis);
    auto dipole_op = sp_basis.quantize(GQCP::Operator::ElectronicDipole());
    dipole_op.transform(rhf.get_C());

    GQCP::Vector<double, 3> total_dipole_moment = GQCP::Operator::NuclearDipole(CO).value() + GQCP::calculateElectronicDipoleMoment(dipole_op, D);
    BOOST_CHECK(std::abs(total_dipole_moment.norm() - (0.049)) < 1.0e-03);
}


BOOST_AUTO_TEST_CASE ( dipole_N2_STO_3G ) {

    // Check that the dipole moment of N2 is zero

    // Initialize the molecule and molecular Hamiltonian parameters for N2
    GQCP::Nucleus N_1 (7, 0.0, 0.0, 0.0);
    GQCP::Nucleus N_2 (7, 0.0, 0.0, GQCP::units::angstrom_to_bohr(1.134));  // from CCCBDB, STO-3G geometry
    std::vector<GQCP::Nucleus> nuclei {N_1, N_2};
    GQCP::Molecule N2 (nuclei);

    auto ao_basis = std::make_shared<GQCP::ScalarBasis<GQCP::GTOShell>>(N2, "STO-3G");
    auto ao_mol_ham_par = GQCP::HamiltonianParameters<double>::Molecular(ao_basis);

    size_t K = ao_basis->numberOfBasisFunctions();
    size_t N = N2.numberOfElectrons();

    // Solve the SCF equations
    GQCP::PlainRHFSCFSolver plain_scf_solver (ao_mol_ham_par, N2);  // The DIIS SCF solver seems to find a wrong minimum, so use a plain solver instead
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    double total_energy = rhf.get_electronic_energy() + GQCP::Operator::NuclearRepulsion(N2).value();
    BOOST_REQUIRE(std::abs(total_energy - (-107.500654)) < 1.0e-05);  // from CCCBDB, require a correct RHF solution to be found


    // Calculate the RHF 1-RDM in MO basis
    auto D = GQCP::calculateRHF1RDM(K, N);
    auto D_AO = GQCP::calculateRHFAO1RDM(rhf.get_C(), N);

    // Calculate the dipole integrals, and transform them to the MO basis
    const GQCP::SingleParticleBasis<double, GQCP::GTOShell> sp_basis (*ao_basis);
    auto dipole_op = sp_basis.quantize(GQCP::Operator::ElectronicDipole());
    dipole_op.transform(rhf.get_C());

    GQCP::Vector<double, 3> total_dipole_moment = GQCP::Operator::NuclearDipole(N2).value() + GQCP::calculateElectronicDipoleMoment(dipole_op, D);
    BOOST_CHECK(std::abs(total_dipole_moment.norm() - (0.0)) < 1.0e-08);
}


/**
 *  Test the Dyson amplitude algorithm against manually calculated amplitudes for two (normalized) toy wave functions
 */
BOOST_AUTO_TEST_CASE ( dyson_amplitudes) {

    GQCP::VectorX<double> reference_amplitudes = GQCP::VectorX<double>::Zero(2); 
    reference_amplitudes << 0.46752445464799997, 0.8415438533449999;

    const size_t K = 2;
    const size_t N = 2;
 
    const GQCP::ProductFockSpace fock_space1 (K, N/2, N/2);
    const GQCP::ProductFockSpace fock_space2 (K, N/2, N/2-1);
    const GQCP::ProductFockSpace fock_space3 (K, N/2-1, N/2);
  
    GQCP::VectorX<double> vec1 = GQCP::VectorX<double>::Zero(4); 
    vec1 << 0.182574, 0.365148, 0.547723, 0.730297;
    GQCP::VectorX<double> vec2 = GQCP::VectorX<double>::Zero(2); 
    vec2 << 0.640184, 0.768221;
    
    const auto wavefunction1 = GQCP::WaveFunction(fock_space1, vec1);
    const auto wavefunction2 = GQCP::WaveFunction(fock_space2, vec2);
    const auto wavefunction3 = GQCP::WaveFunction(fock_space3, vec2);

    const auto dyson_coefficients_beta = GQCP::calculateDysonAmplitudes(wavefunction1, wavefunction2);  // coefficients with a difference in beta occupance
    const auto dyson_coefficients_alpha = GQCP::calculateDysonAmplitudes(wavefunction1, wavefunction3);  // coefficients with a difference in alpha occupance

    BOOST_CHECK(dyson_coefficients_beta.isApprox(reference_amplitudes, 1.0e-6));
    BOOST_CHECK(dyson_coefficients_alpha.isApprox(reference_amplitudes, 1.0e-6));
}
