#include "Basis/ScalarBasis/GTOBasisSet.hpp"
#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Basis/ScalarBasis/ShellSet.hpp"
#include "Basis/SpinorBasis/GSpinorBasis.hpp"
#include "Mathematical/Functions/CartesianGTO.hpp"
#include "Mathematical/Functions/EvaluableLinearCombination.hpp"
#include "Operator/FirstQuantized/CurrentDensityOperator.hpp"
#include "Operator/FirstQuantized/LinearMomentumOperator.hpp"
#include "Operator/FirstQuantized/OverlapOperator.hpp"

#include <iostream>


int main() {

    // Create an STO-3G basisset on (a toy geometry of) H2O.
    const GQCP::Nucleus h1 {1, 0.0, 0.0, 0.0};
    const GQCP::Nucleus o {8, 0.0, 0.0, 1.0};
    const GQCP::Nucleus h2 {1, 0.0, 0.0, 2.0};
    const GQCP::Molecule molecule {{h1, o, h2}};

    GQCP::RSpinOrbitalBasis<GQCP::complex, GQCP::GTOShell> spin_orbital_basis {molecule, "STO-3G"};
    // const GQCP::HomogeneousMagneticField B {{0.0, 0.0, 1.0}};  // Gauge origin at the origin.
    // const GQCP::RSpinOrbitalBasis<GQCP::complex, GQCP::LondonGTOShell> london_spin_orbital_basis {molecule, "STO-3G", B};

    // auto j_op = spin_orbital_basis.quantize(GQCP::CurrentDensityOperator());
    // auto j_op = london_spin_orbital_basis.quantize(GQCP::CurrentDensityOperator());
    auto p_op = spin_orbital_basis.quantize(GQCP::LinearMomentumOperator());
    // auto s_op = spinor_basis.quantize(GQCP::OverlapOperator());

    // auto j = j_op.allParameters()[0];
    auto p = p_op.allParameters()[0];
    // auto s = s_op.parameters();

    // initialize vector
    const GQCP::Vector<double, 3> r = {0.1, 0.2, 0.3};

    // std::cout << j(0,1)(r) << std::endl;
    std::cout << p(0,1)(r) << std::endl;
    // std::cout << s(0, 1)(r) << std::endl;

    // // select one single AO
    // const auto example_shell = shellset.asVector()[4];
    // const auto basis_functions = example_shell.basisFunctions();

    // // initialize vector
    // const GQCP::Vector<double, 3> r = {0.1, 0.2, 0.3};

    // // Evaluate each basis function at r
    // for (size_t i = 0; i < basis_functions.size(); ++i) {
    //     double value = basis_functions[i](r);  // uses operator() from EvaluableLinearCombination
    //     std::cout << "AO " << i << " value at r = " << r.transpose() << " is " << value << std::endl;
    // }

    return 0;
}