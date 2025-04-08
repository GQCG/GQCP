#include "Basis/ScalarBasis/GTOBasisSet.hpp"
#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Basis/ScalarBasis/ShellSet.hpp"
#include "Mathematical/Functions/CartesianGTO.hpp"
#include "Mathematical/Functions/EvaluableLinearCombination.hpp"
#include <iostream>
#include "Operator/FirstQuantized/CurrentDensityOperator.hpp"
#include "Basis/SpinorBasis/RSpinOrbitalBasis.hpp"


int main() {

    // Create an STO-3G basisset on (a toy geometry of) H2O.
    const GQCP::Nucleus h1 {1, 0.0, 0.0, 0.0};
    const GQCP::Nucleus o {8, 0.0, 0.0, 1.0};
    const GQCP::Nucleus h2 {1, 0.0, 0.0, 2.0};
    const GQCP::Molecule molecule {{h1, o, h2}};

    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "STO-3G"};
    auto j_op = spin_orbital_basis.quantize(GQCP::CurrentDensityOperator());

    auto j = j_op.allParameters()[0];

    // initialize vector
    const GQCP::Vector<double, 3> r = {0.1, 0.2, 0.3};

    std::cout << j(0,1)(r) << std::endl;

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