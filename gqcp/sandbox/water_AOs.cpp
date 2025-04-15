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
    const GQCP::Nucleus o {8, 1.0, 0.0, 0.0};
    const GQCP::Nucleus h2 {1, 2.0, 0.0, 0.0};
    const GQCP::Molecule molecule {{h1, o, h2}};

    const auto B = GQCP::HomogeneousMagneticField {{0.0, 0.0, 1.0}};  // Gauge origin at cartesian origin.
    // const auto B = GQCP::HomogeneousMagneticField {{0.0, 0.0, 1.0}, {0.4, 12.2, -0.789}};  // Gauge origin at random point in space.
    auto spin_orbital_basis = GQCP::RSpinOrbitalBasis<GQCP::complex, GQCP::LondonGTOShell> {molecule, "STO-3G", B};
    
    // gives all basis functions in the basis. vector of (contraction coeff, basis function) pairs.
    // these are NOT necessarily the AO basis functions, although by default the AO basis is the expansion basis.
    auto spatial_basis_functions_vector = spin_orbital_basis.spatialOrbitals(); 
    print_type(spatial_basis_functions_vector);

    // identical:
    std::cout << "amount of orbitals (MO level)" << std::endl;
    std::cout << spatial_basis_functions_vector.size() << std::endl; 
    std::cout << spin_orbital_basis.numberOfSpatialOrbitals() << std::endl; 

    // select one single orbital
    auto orbital1 = spatial_basis_functions_vector[0]; 
    std::cout << "type of spin_orbital_basis.spatialOrbitals()[0]: ";
    print_type(orbital1);

    // this orbital is in its turn expanded in the AO basis.
    // eg for H2O STO-3G, this is 2x 1 for H, 3 for O = 5 total.
    const auto AOs = orbital1.functions();
    std::cout << "type of spin_orbital_basis.spatialOrbitals()[0].functions(): ";
    print_type(AOs);

    // the AOs are LAOs in this case
    auto example_LAO = AOs[0];
    std::cout << "type of spin_orbital_basis.spatialOrbitals()[0].functions()[0]: ";
    print_type(example_LAO);

    // these can be evaluated at a point in space:
    GQCP::Vector<double, 3> r = {1.3, -0.9, 3.7};
    std::cout << "example LAO value at" << r.transpose() << ": " << example_LAO(r) << std::endl;

    // this LAO is a contraction of three London GTOs. These can be indiviually evaluated and have a corresponding phase factor etc.
    auto london_primitives = example_LAO.functions();
    std::cout << "the LAOs are made up of " << london_primitives.size() << " primitives" << std::endl;
    std::cout << "example phase factor" << london_primitives[0].phaseFactor(r) << std::endl;
    // from the primitves, we can determine the center of the AO
    std::cout << "example LAO is centered at " << example_LAO.functions()[0].cartesianGTO().center().transpose() << std::endl;

    // let's loop through all AOs, and gather their information.
    const int n_ao = AOs.size();
    for (size_t i = 0; i < n_ao; i++) {
        std::cout << i << std::endl;
        // gather relevant AO
        auto AO_i = AOs[i];
        // TODO: ensure to embedNormalizationFactorsOfPrimitives()

        // evaluate its value at r
        std::cout << "value at r: " << AO_i(r) << std::endl;
        // from primitves, get more info.
        auto primitives = AO_i.functions();
        // phase factor at r
        std::cout << "phase factor at r: " << primitives[0].phaseFactor(r) << std::endl;
        // corresponding nucleus origin?
        std::cout << "origin at " << primitives[0].cartesianGTO().center().transpose() << std::endl;
    }

    // note: see SimpleSpinOrbitalBasis.hpp: the normalization factors of the primitives are already embedded in the contraction coefficients of the underlying shells.

    return 0;
}