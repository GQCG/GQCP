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
#include "Basis/transform.hpp"
#include "Operator/FirstQuantized/Operator.hpp"
#include "Processing/Properties/RHFElectricalResponseSolver.hpp"
#include "Processing/Properties/properties.hpp"
#include "Processing/RDM/RDMCalculator.hpp"
#include "QCMethod/CI/HamiltonianBuilder/FCI.hpp"
#include "QCMethod/HF/RHF/DiagonalRHFFockMatrixObjective.hpp"
#include "QCMethod/HF/RHF/RHF.hpp"
#include "QCMethod/HF/RHF/RHFSCFSolver.hpp"
#include "Utilities/units.hpp"


BOOST_AUTO_TEST_CASE(dipole_CO_STO_3G) {

    // Initialize the molecule and molecular Hamiltonian for CO
    GQCP::Nucleus C {6, 0.0, 0.0, 0.0};
    GQCP::Nucleus O {8, 0.0, 0.0, GQCP::units::angstrom_to_bohr(1.145)};  // from CCCBDB, STO-3G geometry
    GQCP::Molecule CO {{C, O}};

    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis {CO, "STO-3G"};
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, CO);  // in an AO basis

    size_t K = spinor_basis.numberOfSpatialOrbitals();
    size_t N = CO.numberOfElectrons();

    // Solve the SCF equations
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(CO.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters());
    auto diis_rhf_scf_solver = GQCP::RHFSCFSolver<double>::DIIS();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective {sq_hamiltonian};
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


BOOST_AUTO_TEST_CASE(dipole_N2_STO_3G) {

    // Check that the dipole moment of N2 is zero

    // Initialize the molecule and the molecular Hamiltonian for N2
    GQCP::Nucleus N_1 {7, 0.0, 0.0, 0.0};
    GQCP::Nucleus N_2 {7, 0.0, 0.0, GQCP::units::angstrom_to_bohr(1.134)};  // from CCCBDB, STO-3G geometry
    GQCP::Molecule N2 {{N_1, N_2}};

    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis {N2, "STO-3G"};
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, N2);  // in an AO basis

    size_t K = spinor_basis.numberOfSpatialOrbitals();
    size_t N = N2.numberOfElectrons();

    // Solve the SCF equations
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(N2.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective {sq_hamiltonian};
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
BOOST_AUTO_TEST_CASE(h2_polarizability_RHF) {

    // Initialize the reference value
    const double ref_alpha_zz = 1.08428;


    // Initialize the molecule and the Hamiltonian in the AO basis
    GQCP::Nucleus H1 {1, 0.0, 0.0, 0.0};
    GQCP::Nucleus H2 {1, 0.0, 0.0, 0.5};
    GQCP::Molecule h2 {{H1, H2}, 0};

    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis(h2, "STO-3G");
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, h2);  // in the AO basis


    // Do the RHF calculation to get the canonical RHF orbitals
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(h2.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective(sq_hamiltonian);
    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters();


    // Transform the orbitals to the RHF basis and prepare the dipole integrals in the RHF basis
    GQCP::basisTransform(spinor_basis, sq_hamiltonian, rhf_parameters.coefficientMatrix());
    const auto dipole_op = spinor_basis.quantize(GQCP::Operator::ElectronicDipole());


    // Find the RHF wave function response
    GQCP::RHFElectricalResponseSolver cphf_solver(h2.numberOfElectrons() / 2);
    const auto x = cphf_solver.calculateWaveFunctionResponse(sq_hamiltonian, dipole_op);


    // Calculate the RHF polarizability
    const auto F_p = cphf_solver.calculateParameterResponseForce(dipole_op);
    const auto alpha = GQCP::calculateElectricPolarizability(F_p, x);
    const auto alpha_zz = alpha(2, 2);


    BOOST_CHECK(std::abs(alpha_zz - ref_alpha_zz) < 1.0e-05);
}


/**
 *  Test the Dyson algorithm against manually calculated coefficients for two (normalized) toy wave functions
 */
BOOST_AUTO_TEST_CASE(dyson_coefficients) {

    // Set up the manually calculated references
    GQCP::VectorX<double> reference_amplitudes_beta = GQCP::VectorX<double>::Zero(2);
    reference_amplitudes_beta << 0.537653264399, 0.794791398869;

    GQCP::VectorX<double> reference_amplitudes_alpha = GQCP::VectorX<double>::Zero(2);
    reference_amplitudes_alpha << 0.39739531532399996, 0.9116729926689999;

    // Set up the toy wave functions
    const size_t K = 2;
    const size_t N = 2;

    const GQCP::SpinResolvedONVBasis fock_space1 {K, N / 2, N / 2};
    const GQCP::SpinResolvedONVBasis fock_space2 {K, N / 2, N / 2 - 1};
    const GQCP::SpinResolvedONVBasis fock_space3 {K, N / 2 - 1, N / 2};

    GQCP::VectorX<double> coeffs1 = GQCP::VectorX<double>::Zero(4);
    coeffs1 << 0.182574, 0.365148, 0.547723, 0.730297;
    GQCP::VectorX<double> coeffs2 = GQCP::VectorX<double>::Zero(2);
    coeffs2 << 0.640184, 0.768221;

    const auto linear_expansion1 = GQCP::LinearExpansion<GQCP::SpinResolvedONVBasis>(fock_space1, coeffs1);
    const auto linear_expansion2 = GQCP::LinearExpansion<GQCP::SpinResolvedONVBasis>(fock_space2, coeffs2);
    const auto linear_expansion3 = GQCP::LinearExpansion<GQCP::SpinResolvedONVBasis>(fock_space3, coeffs2);

    // Calculate the coefficients of the Dyson orbitals and check with the reference
    const auto dyson_coefficients_beta = GQCP::calculateDysonOrbitalCoefficients(linear_expansion1, linear_expansion2);   // coefficients with a difference in beta occupation
    const auto dyson_coefficients_alpha = GQCP::calculateDysonOrbitalCoefficients(linear_expansion1, linear_expansion3);  // coefficients with a difference in alpha occupation

    BOOST_CHECK(dyson_coefficients_beta.isApprox(reference_amplitudes_beta, 1.0e-6));
    BOOST_CHECK(dyson_coefficients_alpha.isApprox(reference_amplitudes_alpha, 1.0e-6));
}

#include "Basis/Integrals/IntegralCalculator.hpp"
#include "Mathematical/Grid/CubicGrid.hpp"
#include "Mathematical/Optimization/LinearEquation/LinearEquationEnvironment.hpp"
#include "Mathematical/Optimization/LinearEquation/LinearEquationSolver.hpp"
#include "Mathematical/Representation/LeviCivitaTensor.hpp"
#include "Utilities/aliases.hpp"
#include "Utilities/literals.hpp"

const GQCP::LeviCivitaTensor epsilon {};


using namespace GQCP::literals;


BOOST_AUTO_TEST_CASE(current_sandbox) {

    // Set up the molecular Hamiltonian in the AO spin-orbital basis.
    const auto molecule = GQCP::Molecule::HChain(2, 1.0);
    const auto N_P = molecule.numberOfElectronPairs();

    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis {molecule, "6-31G**"};
    const auto K = spinor_basis.numberOfSpatialOrbitals();

    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, molecule);  // in the AO basis


    // Find the RHF canonical orbitals and transform the integrals to that basis.
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(molecule.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters());
    auto rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective(sq_hamiltonian);
    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, rhf_scf_solver, rhf_environment).groundStateParameters();

    GQCP::basisTransform(spinor_basis, sq_hamiltonian, rhf_parameters.coefficientMatrix());


    /*
     *  SET UP THE LINEAR RESPONSE EQUATIONS.
     */
    const auto& scalar_basis = spinor_basis.scalarBasis();
    const auto orbital_space = rhf_parameters.orbitalSpace();

    const auto dim = orbital_space.numberOfExcitations(GQCP::OccupationType::k_occupied, GQCP::OccupationType::k_virtual);
    std::cout << orbital_space.description() << std::endl;


    // Calculate the response force constant for RHF.
    const auto k_kappa = rhf_parameters.calculateOrbitalHessianForImaginaryResponse(sq_hamiltonian);

    std::cout << "k_kappa: " << std::endl
              << k_kappa.asMatrix() << std::endl
              << std::endl;


    // Calculate the response force for the magnetic field perturbation.
    auto angular_momentum_engine = GQCP::IntegralEngine::InHouse(GQCP::Operator::AngularMomentum());  // gauge origin == zero
    const auto L = GQCP::IntegralCalculator::calculate(angular_momentum_engine, scalar_basis.shellSet(), scalar_basis.shellSet());

    // TODO: ! TRANSFORM TO CANONICAL ORBITAL BASIS ! //

    GQCP::Matrix<GQCP::complex, GQCP::Dynamic, 3> F_kappa_B = GQCP::Matrix<GQCP::complex, GQCP::Dynamic, 3>::Zero(dim, 3);

    for (size_t m = 0; m < 3; m++) {
        auto F_kappa_B_m = orbital_space.initializeRepresentableObjectFor<GQCP::complex>(GQCP::OccupationType::k_virtual, GQCP::OccupationType::k_occupied);

        for (const auto& a : orbital_space.indices(GQCP::OccupationType::k_virtual)) {
            for (const auto& i : orbital_space.indices(GQCP::OccupationType::k_occupied)) {
                F_kappa_B_m(a, i) = -2.0 * L[m](i, a);
            }
        }
        F_kappa_B.col(m) = F_kappa_B_m.asVector();
    }


    std::cout << "F_kappa_B: " << std::endl
              << F_kappa_B << std::endl
              << std::endl;


    // Calculate the response force for the gauge origin perturbation.
    const auto p = spinor_basis.quantize(GQCP::Operator::LinearMomentum()).allParameters();  // the linear momentum integrals in the canonical RHF basis


    GQCP::Matrix<GQCP::complex, GQCP::Dynamic, 6> F_kappa_G = GQCP::Matrix<GQCP::complex, GQCP::Dynamic, 6>::Zero(dim, 6);  // don't include the diagonal xx, yy, zz

    size_t column_index = 0;
    for (size_t m = 0; m < 3; m++) {
        for (size_t n = 0; n < 3; n++) {
            if (m == n) {  // skip the diagonal xx, yy, zz
                continue;
            }

            auto F_kappa_G_mn = orbital_space.initializeRepresentableObjectFor<GQCP::complex>(GQCP::OccupationType::k_virtual, GQCP::OccupationType::k_occupied);

            for (const auto& a : orbital_space.indices(GQCP::OccupationType::k_virtual)) {
                for (const auto& i : orbital_space.indices(GQCP::OccupationType::k_occupied)) {
                    const auto f = epsilon.nonZeroIndex(m, n);

                    F_kappa_G_mn(a, i) = -2.0 * epsilon(m, n, f) * p[f](i, a);
                }
            }

            F_kappa_G.col(column_index) = F_kappa_G_mn.asVector();
            column_index++;
        }
    }


    std::cout << "F_kappa_G: " << std::endl
              << F_kappa_G << std::endl
              << std::endl;


    /*
     *  CALCULATE THE RESPONSES BY SOLVING THE LINEAR RESPONSE EQUATIONS.
     */
    auto environment_B = GQCP::LinearEquationEnvironment<GQCP::complex>(k_kappa.asMatrix(), -F_kappa_B);
    auto solver_B = GQCP::LinearEquationSolver<GQCP::complex>::HouseholderQR();
    solver_B.perform(environment_B);

    const auto x_B = environment_B.x;

    auto environment_G = GQCP::LinearEquationEnvironment<GQCP::complex>(k_kappa.asMatrix(), -F_kappa_G);
    auto solver_G = GQCP::LinearEquationSolver<GQCP::complex>::HouseholderQR();
    solver_G.perform(environment_G);

    const auto x_G = environment_G.x;


    std::cout << "x_B: " << std::endl
              << x_B << std::endl
              << std::endl;
    std::cout << "x_G: " << std::endl
              << x_G << std::endl
              << std::endl;


    /*
     *  EVALUATE THE FIELD-FREE CURRENT DENSITY MATRIX ELEMENTS
     */
    // We require the spatial orbitals and their gradients.
    const auto spatial_orbitals = spinor_basis.spatialOrbitals();

    using Primitive = GQCP::CartesianGTO;                                   // the primitives are Cartesian GTOs
    using BasisFunction = GQCP::LinearCombination<double, Primitive>;       // the basis functions are contracted GTOs
    using SpatialOrbital = GQCP::LinearCombination<double, BasisFunction>;  // spatial orbitals are linear combinations of basis functions


    using PrimitiveDerivative = GQCP::LinearCombination<double, Primitive>;  // the derivative of a Cartesian GTO is a linear combination of Cartesian GTOs
    using BasisFunctionDerivative = GQCP::LinearCombination<double, PrimitiveDerivative>;
    using SpatialOrbitalDerivative = GQCP::LinearCombination<double, BasisFunctionDerivative>;

    // Calculate the gradients of the spatial orbitals.
    // TODO: fix RSpinorBasis_test duration: check what happened in LinearCombination
    std::vector<GQCP::Vector<SpatialOrbitalDerivative, 3>> spatial_orbital_gradients {K};
    for (size_t m = 0; m < 3; m++) {
        for (size_t p = 0; p < K; p++) {
            std::cout << "p: " << p << std::endl;
            const auto& spatial_orbital = spatial_orbitals[p];

            // A spatial orbital is a linear combination of basis functions (which are contracted GTOs).
            const auto& expansion_coefficients = spatial_orbital.coefficients();
            const auto& basis_functions = spatial_orbital.functions();

            std::vector<GQCP::Vector<BasisFunctionDerivative, 3>> basis_function_gradients {K};
            for (size_t mu = 0; mu < K; mu++) {
                std::cout << "mu: " << mu << std::endl;
                const auto& expansion_coefficient = expansion_coefficients[mu];
                const auto& basis_function = basis_functions[mu];
                std::cout << "basis function: " << basis_function.description() << std::endl;
                const auto contraction_length = basis_function.length();

                // A basis function (a.k.a a contracted GTO) is a contraction of Cartesian GTOs.
                const auto& contraction_coefficients = basis_function.coefficients();
                const auto& primitives = basis_function.functions();

                for (size_t d = 0; d < contraction_length; d++) {
                    const auto& contraction_coefficient = contraction_coefficients[d];
                    const auto primitive_gradient = primitives[d].calculateGradient();

                    std::cout << "d: " << d << std::endl;

                    basis_function_gradients[mu](m).append(contraction_coefficient, primitive_gradient(m));
                }

                spatial_orbital_gradients[p](m).append(expansion_coefficient, basis_function_gradients[mu](m));
            }
        }
    }

    std::cout << spatial_orbital_gradients[0](2).description() << std::endl
              << std::endl;


    /*
     *  EVALUATE THE MAGNETIC INDUCIBILITY USING THE IPSOCENTRIC CSGT.
     */
    const auto grid = GQCP::CubicGrid::Centered(GQCP::Vector<double, 3>::Zero(), 10, 0.1);
    std::vector<GQCP::Matrix<GQCP::complex, 3, 3>> J_field_values;
    J_field_values.reserve(1000);

    grid.forEach([&orbital_space, &x_B, &x_G, &J_field_values, &spatial_orbitals, &spatial_orbital_gradients](const GQCP::Vector<double, 3>& r) {
        for (size_t u = 0; u < 3; u++) {
            for (size_t m = 0; m < 3; m++) {
                GQCP::Matrix<GQCP::complex, 3, 3> J = GQCP::Matrix<GQCP::complex, 3, 3>::Zero();

                auto D_m = orbital_space.initializeRepresentableObjectFor<GQCP::complex>(GQCP::OccupationType::k_virtual, GQCP::OccupationType::k_occupied);
                auto x_m_matrix = GQCP::MatrixX<GQCP::complex>::FromColumnMajorVector(x_B.col(m), orbital_space.numberOfOrbitals(GQCP::OccupationType::k_virtual), orbital_space.numberOfOrbitals(GQCP::OccupationType::k_occupied));
                auto x_m = orbital_space.createRepresentableObjectFor<GQCP::complex>(GQCP::OccupationType::k_virtual, GQCP::OccupationType::k_occupied, x_m_matrix);

                for (const auto& a : orbital_space.indices(GQCP::OccupationType::k_virtual)) {
                    for (const auto& i : orbital_space.indices(GQCP::OccupationType::k_occupied)) {
                        GQCP::complex left_value;
                        GQCP::complex right_value;

                        left_value = spatial_orbitals[i](r) * spatial_orbital_gradients[a](u)(r) - spatial_orbital_gradients[i](u)(r) * spatial_orbitals[a](r);

                        right_value = -1.0_ii * x_m(a, i);
                        for (size_t n = 0; n != m; n++) {
                            const auto row_major_index = 3 * m + n;
                            const auto mn = row_major_index < 4 ? row_major_index - 1 : row_major_index - 2;  // The compount index 'mn'.

                            auto epsilon_mn_matrix = GQCP::MatrixX<GQCP::complex>::FromColumnMajorVector(x_G.col(mn), orbital_space.numberOfOrbitals(GQCP::OccupationType::k_virtual), orbital_space.numberOfOrbitals(GQCP::OccupationType::k_occupied));
                            auto epsilon_mn = orbital_space.createRepresentableObjectFor<GQCP::complex>(GQCP::OccupationType::k_virtual, GQCP::OccupationType::k_occupied, epsilon_mn_matrix);

                            right_value += -1.0_ii * epsilon_mn(a, i) * r(n);  // d = r - G_0, and we have used G_0 = O; this is the ipsocentric CSGT step
                        }

                        J(u, m) += -2_ii * left_value * right_value;
                        assert(J(u, m).imag() < 1.0e-12);
                    }
                }
                J_field_values.push_back(J);
            }
        }
    });

    std::cout << "J_field_values: " << std::endl;
    for (const auto& each : J_field_values) {
        std::cout << each << ' ';
    }
    std::cout << std::endl
              << std::endl;


    GQCP::Field<GQCP::Matrix<GQCP::complex, 3, 3>> J {J_field_values};


    /*
     *  CALCULATE THE NICSD_ZZ VALUES.
     */
    const double prefactor = std::pow(1.0 / 137.0, 2);
    std::vector<GQCP::complex> NICSD_zz_values;
    NICSD_zz_values.reserve(1000);

    const auto positions = grid.points();
    const GQCP::Vector<double, 3> R_K = {0.003, 0.003, 0.003};  // reference point for NICSD_zz(R_K)

    for (size_t i = 0; i < positions.size(); i++) {
        const auto& r = positions[i];

        GQCP::complex value {};
        value += J_field_values[i](0, 2) * (r - R_K)(1) - J_field_values[i](1, 2) * (r - R_K)(0);

        NICSD_zz_values.push_back(prefactor * value / std::pow((r - R_K).norm(), 3));
    }


    std::cout << "NICSD_zz_values: " << std::endl;
    for (const auto& each : NICSD_zz_values) {
        std::cout << each << ' ';
    }
    std::cout << std::endl
              << std::endl;


    GQCP::Field<GQCP::complex> NICSD_zz {NICSD_zz_values};
    const auto NICS_zz = grid.integrate(NICSD_zz);

    std::cout << "NICS_zz (at {0.003, 0.003, 0.003}): " << NICS_zz << std::endl;
    // NICS_zz (at {0.003, 0.003, 0.003}): (-1.26902e-08,9.92093e-26)
}
