#include "Basis/SpinorBasis/RSpinOrbitalBasis.hpp"
#include "QCMethod/HF/RHF/DiagonalRHFFockMatrixObjective.hpp"
#include "QCMethod/HF/RHF/RHF.hpp"
#include "QCMethod/HF/RHF/RHFSCFSolver.hpp"
#include "QCModel/HF/RHF.hpp"

// from RHF_test

int main() {

    using namespace GQCP::literals;

    // Set up the molecular Hamiltonian in AO basis.
    const GQCP::Molecule molecule {{GQCP::Nucleus(1, 0.0, 0.0, 0.0), GQCP::Nucleus(1, 0.0, 0.0, 1.0)}};
    const auto N_P = molecule.numberOfElectronPairs();

    const std::string basis_set {"STO-3G"};
    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, basis_set};

    auto hamiltonian = spin_orbital_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));


    // Read in the GAMESS-UK RHF wave function model parameters. Even though GQCP finds the same orbitals and orbital energies, the phase factors (and orbital coefficients for degenerated orbitals) are unlikely to be reproducible.
    GQCP::MatrixX<double> C_matrix {2, 2};
    // clang-format off
    C_matrix << 0.5275464665,  1.5678230259,
                0.5275464665, -1.5678230259;
    // clang-format on
    GQCP::RTransformation<double> C {C_matrix};

    GQCP::VectorX<double> orbital_energies {2};
    orbital_energies << -0.6757801904, 0.9418115528;

    const GQCP::QCModel::RHF<double> rhf_parameters {N_P, orbital_energies, C};


    // Since we're going to work with complex operators, we have to let a complex spin-orbital basis do the quantization.
    GQCP::RSpinOrbitalBasis<GQCP::complex, GQCP::GTOShell> complex_spin_orbital_basis {molecule, basis_set};
    GQCP::RTransformation<GQCP::complex> C_complex {C_matrix.cast<GQCP::complex>()};
    complex_spin_orbital_basis.transform(C_complex);

    spin_orbital_basis.transform(C);
    hamiltonian.transform(C);


    // Calculate the orbital Hessian.
    const auto orbital_space = rhf_parameters.orbitalSpace();
    auto A = rhf_parameters.calculateOrbitalHessianForImaginaryResponse(hamiltonian, orbital_space);


    // Solve the CPHF equations for the angular momentum operator.
    const auto L = complex_spin_orbital_basis.quantize(GQCP::AngularMomentumOperator());
    const auto F_B = rhf_parameters.calculateMagneticFieldResponseForce(L);


    auto environment_B = GQCP::LinearEquationEnvironment<GQCP::complex>(A.asMatrix(), -F_B);
    auto solver_B = GQCP::LinearEquationSolver<GQCP::complex>::HouseholderQR();
    solver_B.perform(environment_B);
    const auto x = environment_B.x;


    // Solve the CPHF equations for the linear momentum operator.
    const auto p = complex_spin_orbital_basis.quantize(GQCP::LinearMomentumOperator());
    const auto F_G = rhf_parameters.calculateGaugeOriginTranslationResponseForce(p);

    // In order to check with the reference values, we have to convert our dyadic Cartesian (i.e. xy, xz, etc.) representation to a Cartesian (i.e. x,y,z) one.
    GQCP::MatrixX<GQCP::complex> F_G_reduced {1, 3};
    F_G_reduced.col(0) = 2 * F_G.col(3);  // x <--> yz
    F_G_reduced.col(1) = 2 * F_G.col(4);  // y <--> zx
    F_G_reduced.col(2) = 2 * F_G.col(0);  // z <--> xy


    auto environment_G = GQCP::LinearEquationEnvironment<GQCP::complex>(A.asMatrix(), -F_G);
    auto solver_G = GQCP::LinearEquationSolver<GQCP::complex>::HouseholderQR();
    solver_G.perform(environment_G);

    const auto y = environment_G.x;

    // In order to check with the reference values, we have to convert our dyadic Cartesian (i.e. xy, xz, etc.) representation to a Cartesian (i.e. x, y, z) one.
    GQCP::MatrixX<GQCP::complex> y_reduced {1, 3};
    y_reduced.col(0) = 2 * y.col(3);  // x <--> yz
    y_reduced.col(1) = 2 * y.col(4);  // y <--> zx
    y_reduced.col(2) = 2 * y.col(0);  // z <--> xy


    // Calculate the ipsocentric magnetic inducibility on a grid and check the results.
    const auto j_op = complex_spin_orbital_basis.quantize(GQCP::CurrentDensityOperator());

    const GQCP::Vector<double, 3> origin {-2.0, -2.0, -2.0};
    const std::array<size_t, 3> steps {6, 6, 6};
    const std::array<double, 3> step_sizes {0.8, 0.8, 0.8};
    const GQCP::CubicGrid grid {origin, steps, step_sizes};

    const auto J_field = GQCP::QCModel::RHF<GQCP::complex>::calculateIpsocentricMagneticInducibility(grid, orbital_space, x, y, j_op);
    const auto& J_values = J_field.values();

    return 0;
}