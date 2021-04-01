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

#define BOOST_TEST_MODULE "properties"

#include <boost/test/unit_test.hpp>

#include "Basis/Integrals/Interfaces/LibintInterfacer.hpp"
#include "Basis/Transformations/transform.hpp"
#include "Operator/FirstQuantized/NuclearDipoleOperator.hpp"
#include "Operator/FirstQuantized/NuclearRepulsionOperator.hpp"
#include "Physical/units.hpp"
#include "Processing/Properties/RHFElectricalResponseSolver.hpp"
#include "Processing/Properties/properties.hpp"
#include "QCMethod/HF/RHF/DiagonalRHFFockMatrixObjective.hpp"
#include "QCMethod/HF/RHF/RHF.hpp"
#include "QCMethod/HF/RHF/RHFSCFSolver.hpp"
#include "Utilities/complex.hpp"


/**
 *  Check the calculation of the CO dipole moment from a CCCBDB reference value.
 */
BOOST_AUTO_TEST_CASE(dipole_CO_STO_3G) {

    // Initialize the molecule and molecular Hamiltonian for CO.
    const GQCP::Nucleus C {6, 0.0, 0.0, 0.0};
    const GQCP::Nucleus O {8, 0.0, 0.0, GQCP::units::angstrom_to_bohr(1.145)};  // From CCCBDB, STO-3G geometry.
    const GQCP::Molecule molecule {{C, O}};

    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "STO-3G"};
    auto sq_hamiltonian = spin_orbital_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));  // In the AO basis.

    const auto K = spin_orbital_basis.numberOfSpatialOrbitals();
    const size_t N = molecule.numberOfElectrons();

    // Solve the RHF SCF equations.
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(molecule.numberOfElectrons(), sq_hamiltonian, spin_orbital_basis.overlap().parameters());
    auto diis_rhf_scf_solver = GQCP::RHFSCFSolver<double>::DIIS();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective {sq_hamiltonian};
    const auto rhf_qc_structure = GQCP::QCMethod::RHF<double>().optimize(objective, diis_rhf_scf_solver, rhf_environment);
    const auto rhf_parameters = rhf_qc_structure.groundStateParameters();

    const double total_energy = rhf_qc_structure.groundStateEnergy() + GQCP::NuclearRepulsionOperator(molecule.nuclearFramework()).value();
    BOOST_REQUIRE(std::abs(total_energy - (-111.225)) < 1.0e-02);  // From CCCBDB, require a correct RHF solution to be found.


    // Calculate the RHF 1-DM and the dipole operator in RHF MO basis.
    const auto D = GQCP::QCModel::RHF<double>::calculateOrthonormalBasis1DM(K, N);
    auto dipole_op = spin_orbital_basis.quantize(GQCP::ElectronicDipoleOperator());
    dipole_op.transform(rhf_parameters.expansion());

    // Calculate the RHF total dipole moment in the MO basis and check with the reference value.
    GQCP::Vector<double, 3> total_dipole_moment = GQCP::NuclearDipoleOperator(molecule.nuclearFramework()).value() + dipole_op.calculateExpectationValue(D).asVector();
    BOOST_CHECK(std::abs(total_dipole_moment.norm() - (0.049)) < 1.0e-03);
}


/**
 *  Check the the RHF dipole moment for N2 is zero.
 */
BOOST_AUTO_TEST_CASE(dipole_N2_STO_3G) {

    // Initialize the molecule and the molecular Hamiltonian.
    const GQCP::Nucleus N_1 {7, 0.0, 0.0, 0.0};
    const GQCP::Nucleus N_2 {7, 0.0, 0.0, GQCP::units::angstrom_to_bohr(1.134)};  // From CCCBDB, STO-3G geometry.
    const GQCP::Molecule molecule {{N_1, N_2}};

    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "STO-3G"};
    auto sq_hamiltonian = spin_orbital_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));  // In the AO basis.

    const auto K = spin_orbital_basis.numberOfSpatialOrbitals();
    const auto N = molecule.numberOfElectrons();

    // Solve the RHF SCF equations.
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(molecule.numberOfElectrons(), sq_hamiltonian, spin_orbital_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective {sq_hamiltonian};
    const auto rhf_qc_structure = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment);
    const auto rhf_parameters = rhf_qc_structure.groundStateParameters();

    const double total_energy = rhf_qc_structure.groundStateEnergy() + GQCP::NuclearRepulsionOperator(molecule.nuclearFramework()).value();
    BOOST_REQUIRE(std::abs(total_energy - (-107.500654)) < 1.0e-05);  // From CCCBDB, require a correct RHF solution to be found.


    // Calculate the RHF 1-DM and the dipole operator in RHF MO basis.
    const auto D = GQCP::QCModel::RHF<double>::calculateOrthonormalBasis1DM(K, N);
    auto dipole_op = spin_orbital_basis.quantize(GQCP::ElectronicDipoleOperator());
    dipole_op.transform(rhf_parameters.expansion());

    // Calculate the RHF total dipole moment in the MO basis and check with the reference value.
    GQCP::Vector<double, 3> total_dipole_moment = GQCP::NuclearDipoleOperator(molecule.nuclearFramework()).value() + dipole_op.calculateExpectationValue(D).asVector();
    BOOST_CHECK(std::abs(total_dipole_moment.norm() - (0.0)) < 1.0e-08);
}


/**
 *  Check the calculation of the zz-component of the polarizability for H2 with a reference value from Psi4-numpy.
 * 
 *  Note that the reference value is generated from Psi4-numpy, with a fix for the Fockian matrix.
 */
BOOST_AUTO_TEST_CASE(h2_polarizability_RHF) {

    // Initialize the reference value.
    const double ref_alpha_zz = 1.08428;


    // Initialize the molecule and the Hamiltonian in the AO basis.
    const GQCP::Nucleus H1 {1, 0.0, 0.0, 0.0};
    const GQCP::Nucleus H2 {1, 0.0, 0.0, 0.5};
    const GQCP::Molecule molecule {{H1, H2}, 0};

    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "STO-3G"};
    auto sq_hamiltonian = spin_orbital_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));  // In the AO basis.


    // Do the RHF calculation to get the canonical RHF orbitals.
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(molecule.numberOfElectrons(), sq_hamiltonian, spin_orbital_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective(sq_hamiltonian);
    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters();


    // Transform the orbitals to the RHF basis and prepare the dipole integrals in the RHF basis.
    GQCP::transform(rhf_parameters.expansion(), spin_orbital_basis, sq_hamiltonian);
    const auto dipole_op = spin_orbital_basis.quantize(GQCP::ElectronicDipoleOperator());


    // Find the RHF wave function response.
    GQCP::RHFElectricalResponseSolver cphf_solver {molecule.numberOfElectrons() / 2};
    const auto x = cphf_solver.calculateWaveFunctionResponse(sq_hamiltonian, dipole_op);


    // Calculate the RHF polarizability and check with the reference value.
    const auto F_p = cphf_solver.calculateParameterResponseForce(dipole_op);
    const auto alpha = GQCP::calculateElectricPolarizability(F_p, x);
    const auto alpha_zz = alpha(2, 2);

    BOOST_CHECK(std::abs(alpha_zz - ref_alpha_zz) < 1.0e-05);
}


#include "Mathematical/Grid/CubicGrid.hpp"
#include "Operator/FirstQuantized/NuclearRepulsionOperator.hpp"

#include <iomanip>


BOOST_AUTO_TEST_CASE(benzene_inducibility) {

    std::cout << std::setprecision(15);

    // Set up the molecular Hamiltonian in the AO spin-orbital basis.
    const auto molecule = GQCP::Molecule::HChain(2, 1.0);
    std::cout << molecule.description() << std::endl;
    // const auto molecule = GQCP::Molecule::ReadXYZ("data/benzene.xyz");

    // Set up a cubic grid.
    const GQCP::Vector<double, 3> origin {-2.0, -2.0, -2.0};
    const std::array<size_t, 3> steps {5, 5, 5};
    const std::array<double, 3> step_sizes {0.8, 0.8, 0.8};
    const GQCP::CubicGrid grid {origin, steps, step_sizes};


    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "STO-3G"};
    auto hamiltonian = spin_orbital_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));  // In the AO basis.
    std::cout << "Number of orbitals: " << spin_orbital_basis.numberOfSpatialOrbitals() << std::endl;
    std::cout << "Number of electrons: " << molecule.numberOfElectrons() << std::endl;

    // Do the RHF calculation to get the canonical RHF orbitals.
    auto environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(molecule.numberOfElectrons(), hamiltonian, spin_orbital_basis.overlap().parameters());
    auto solver = GQCP::RHFSCFSolver<double>::DIIS(1.0e-06);
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective {hamiltonian};
    const auto qc_structure = GQCP::QCMethod::RHF<double>().optimize(objective, solver, environment);
    const auto rhf_parameters = qc_structure.groundStateParameters();

    const auto& C = rhf_parameters.expansion();
    std::cout << "RHF orbitals coefficient matrix:" << std::endl
              << C.matrix() << std::endl
              << std::endl;
    // std::cout << "RHF SCF done" << std::endl;
    std::cout << "Total RHF energy: " << GQCP::NuclearRepulsionOperator(molecule.nuclearFramework()).value() + qc_structure.groundStateEnergy() << std::endl;


    // Set up the linear response equations (in the RHF MO basis).
    const GQCP::RTransformation<GQCP::complex> C_complex {C.matrix().cast<GQCP::complex>()};
    const GQCP::RSpinOrbitalBasis<GQCP::complex, GQCP::GTOShell> complex_spin_orbital_basis {spin_orbital_basis.scalarBasis(), C_complex};

    spin_orbital_basis.transform(C);
    hamiltonian.transform(C);
    const auto orbital_space = rhf_parameters.orbitalSpace();

    const auto k_kappa = rhf_parameters.calculateOrbitalHessianForImaginaryResponse(hamiltonian, orbital_space).asMatrix();
    std::cout << "k_kappa (Response force constant matrix A in Ax=-b): " << std::endl
              << k_kappa << std::endl
              << std::endl;
    // std::cout << "k_kappa done" << std::endl;

    const auto L = complex_spin_orbital_basis.quantize(GQCP::AngularMomentumOperator());
    std::cout << "Integrals over the angular momentum operator (x,y,z):" << std::endl;
    for (const auto& l : L.allParameters()) {
        std::cout << l << std::endl
                  << std::endl;
    }
    const auto F_kappa_B = rhf_parameters.calculateMagneticFieldResponseForce(L);
    std::cout << "F_kappa_B (Response force vector b in Ax=-b) for magnetic field perturbation:" << std::endl
              << F_kappa_B << std::endl
              << std::endl;
    // std::cout << "F_kappa_B done" << std::endl;

    const auto p = complex_spin_orbital_basis.quantize(GQCP::LinearMomentumOperator());
    std::cout << "Integrals over the linear momentum operator (x,y,z):" << std::endl;
    for (const auto& P : p.allParameters()) {
        std::cout << P << std::endl
                  << std::endl;
    }
    const auto F_kappa_G = rhf_parameters.calculateGaugeOriginTranslationResponseForce(p);
    std::cout << "F_kappa_G (Response force vector b in Ax=-b) for gauge origin translation perturbation:" << std::endl
              << F_kappa_G << std::endl
              << std::endl;
    // std::cout << "F_kappa_G done" << std::endl;


    // Solve the linear response equations.
    auto environment_B = GQCP::LinearEquationEnvironment<GQCP::complex>(k_kappa, -F_kappa_B);
    auto solver_B = GQCP::LinearEquationSolver<GQCP::complex>::HouseholderQR();
    solver_B.perform(environment_B);

    const auto x_B = environment_B.x;
    std::cout << "x_B (Linear response x in Ax=-b) for magnetic field perturbation:" << std::endl
              << x_B << std::endl
              << std::endl;
    // std::cout << "x_B done" << std::endl;

    auto environment_G = GQCP::LinearEquationEnvironment<GQCP::complex>(k_kappa, -F_kappa_G);
    auto solver_G = GQCP::LinearEquationSolver<GQCP::complex>::HouseholderQR();
    solver_G.perform(environment_G);

    const auto x_G = environment_G.x;
    std::cout << "x_G (Linear response x in Ax=-b) for gauge origin translation perturbation:" << std::endl
              << x_G << std::endl
              << std::endl;
    // std::cout << "x_G done" << std::endl;


    // Evaluate the ipsocentric CSGT magnetic inducibility.
    std::cout << "j_op evaluations" << std::endl;
    const auto j_op = complex_spin_orbital_basis.quantize(GQCP::CurrentDensityOperator());
    // grid.forEach([&j_op](const GQCP::Vector<double, 3>& r) {
    //     std::cout << "r:" << std::endl
    //               << r << std::endl;
    //     std::cout << j_op.evaluate(r).parameters() << std::endl;
    // });

    // std::cout << "j_op done" << std::endl;
    const auto J_field = GQCP::QCModel::RHF<GQCP::complex>::calculateIpsocentricMagneticInducibility(grid, orbital_space, x_B, x_G, j_op);
    // std::cout << "J_field done" << std::endl;

    std::cout << std::endl
              << std::endl;
    const auto J_field_values = J_field.values();


    std::ofstream file {"h2_STO-3G.data"};
    const auto points = grid.points();
    for (size_t index = 0; index < points.size(); index++) {
        file << "Grid point:" << std::endl
             << points[index] << std::endl
             << std::endl;

        file << "Inducibility value (xx,xy,xz,yx,yy,yz,zx,zy,zz):" << std::endl
             << J_field_values[index].real() << std::endl
             << std::endl;
    }
}
