#include "Basis/ScalarBasis/GTOBasisSet.hpp"
#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Basis/ScalarBasis/ShellSet.hpp"
#include "Mathematical/Functions/CartesianGTO.hpp"
#include "Mathematical/Functions/EvaluableLinearCombination.hpp"
#include "Basis/SpinorBasis/RSpinOrbitalBasis.hpp"
#include <iostream>

// for determining type of object
#include <typeinfo>
#include <cxxabi.h>


// Add this helper to demangle and print type info
template<typename T>
void print_type(const T& obj) {
    int status;
    const char* name = typeid(obj).name();
    char* realname = abi::__cxa_demangle(name, nullptr, nullptr, &status);
    std::cout << "Type: " << (status == 0 ? realname : name) << std::endl;
    free(realname);
}


int main() {

    // Create an STO-3G basisset on (a toy geometry of) H2O.
    const GQCP::Nucleus h1 {1, 0.0, 0.0, 0.0};
    const GQCP::Nucleus o {8, 0.0, 0.0, 1.0};
    const GQCP::Nucleus h2 {1, 0.0, 0.0, 2.0};
    const GQCP::Molecule molecule {{h1, o, h2}};

    const auto B = GQCP::HomogeneousMagneticField {{0.0, 0.0, 1.0}};  // Gauge origin at the origin.
    auto spin_orbital_basis = GQCP::RSpinOrbitalBasis<GQCP::complex, GQCP::LondonGTOShell> {molecule, "STO-3G", B};
    print_type(spin_orbital_basis);

    print_type(spin_orbital_basis.scalarBasis().shellSet());

    // spin_orbital_basis.scalarBasis().shellSet()[0].embedNormalizationFactorsOfPrimitives();
    
    // gives all AOs (all shells) in the basis. vector of (contraction coeff, function) pairs
    // eg for H2O STO-3G, this is 2x 1 for H, 3 for O = 5 total.
    auto all_AOs_vector = spin_orbital_basis.spatialOrbitals(); 
    print_type(all_AOs_vector);

    // identical:
    std::cout << all_AOs_vector.size() << std::endl; 
    std::cout << spin_orbital_basis.numberOfSpatialOrbitals() << std::endl; 

    // select one single AO
    auto AO1 = all_AOs_vector[0]; 
    print_type(AO1);

    // example_shell.basisFunctions();
    // example_shell.embedNormalizationFactorsOfPrimitives();

    auto basis_functions = AO1.functions();
    print_type(basis_functions);

    // initialize vector
    GQCP::Vector<double, 3> r = {0.1, 0.2, 0.3};

    // try to call AO value directly
    std::cout << AO1(r) << std::endl;
    std::cout << all_AOs_vector[1](r) << std::endl;

    // Evaluate each basis function at r
    for (size_t i = 0; i < basis_functions.size(); ++i) {
        print_type(basis_functions[i]);
        GQCP::complex value = basis_functions[i](r);  // uses operator() from EvaluableLinearCombination
        std::cout << "AO " << i << " value at r = " << r.transpose() << " is " << value << std::endl;
    }

    return 0;
}