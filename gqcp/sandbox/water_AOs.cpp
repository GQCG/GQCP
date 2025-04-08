#include "Basis/ScalarBasis/GTOBasisSet.hpp"
#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Basis/ScalarBasis/ShellSet.hpp"
#include "Mathematical/Functions/CartesianGTO.hpp"
#include "Mathematical/Functions/EvaluableLinearCombination.hpp"
#include <iostream>


int main() {

    // Create an STO-3G basisset on (a toy geometry of) H2O.
    const GQCP::Nucleus h1 {1, 0.0, 0.0, 0.0};
    const GQCP::Nucleus o {8, 0.0, 0.0, 1.0};
    const GQCP::Nucleus h2 {1, 0.0, 0.0, 2.0};
    const GQCP::Molecule molecule {{h1, o, h2}};

    const auto shellset = GQCP::GTOBasisSet("STO-3G").generate(molecule);

    // select one single AO
    const auto example_shell = shellset.asVector()[4];
    const auto basis_functions = example_shell.basisFunctions();

    // initialize vector
    const GQCP::Vector<double, 3> r = {0.1, 0.2, 0.3};

    // Evaluate each basis function at r
    for (size_t i = 0; i < basis_functions.size(); ++i) {
        double value = basis_functions[i](r);  // uses operator() from EvaluableLinearCombination
        std::cout << "AO " << i << " value at r = " << r.transpose() << " is " << value << std::endl;
    }

    return 0;
}