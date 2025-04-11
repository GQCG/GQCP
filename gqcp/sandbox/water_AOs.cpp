#include "Basis/ScalarBasis/GTOBasisSet.hpp"
#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Basis/ScalarBasis/ShellSet.hpp"
#include "Mathematical/Functions/CartesianGTO.hpp"
#include "Mathematical/Functions/EvaluableLinearCombination.hpp"
#include "Basis/SpinorBasis/RSpinOrbitalBasis.hpp"
#include <iostream>


int main() {

    // Create an STO-3G basisset on (a toy geometry of) H2O.
    const GQCP::Nucleus h1 {1, 0.0, 0.0, 0.0};
    const GQCP::Nucleus o {8, 0.0, 0.0, 1.0};
    const GQCP::Nucleus h2 {1, 0.0, 0.0, 2.0};
    const GQCP::Molecule molecule {{h1, o, h2}};

    const auto B = GQCP::HomogeneousMagneticField {{0.0, 0.0, 1.0}};  // Gauge origin at the origin.
    auto spin_orbital_basis = GQCP::RSpinOrbitalBasis<GQCP::complex, GQCP::LondonGTOShell> {molecule, "STO-3G", B}; 
    
    // gives all AOs (all shells) in the basis. vector of (contraction coeff, function) pairs
    // eg for H2O STO-3G, this is 2x 1 for H, 3 for O = 5 total.
    auto shellset = spin_orbital_basis.spatialOrbitals(); 

    // select one single AO
    auto example_shell = shellset[0]; 

    std::cout << shellset.size() << std::endl; 

    std::cout << spin_orbital_basis.numberOfSpatialOrbitals() << std::endl; 

    example_shell.basisFunctions();
    // example_shell.embedNormalizationFactorsOfPrimitives();

    auto basis_functions = example_shell.functions();

    // initialize vector
    GQCP::Vector<double, 3> r = {0.1, 0.2, 0.3};

    // Evaluate each basis function at r
    for (size_t i = 0; i < basis_functions.size(); ++i) {
        GQCP::complex value = basis_functions[i](r);  // uses operator() from EvaluableLinearCombination
        std::cout << "AO " << i << " value at r = " << r.transpose() << " is " << value << std::endl;
    }

    return 0;
}