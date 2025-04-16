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


int main() {

    // Create an STO-3G basisset on (a toy geometry of) H2O.
    const GQCP::Nucleus h1 {1, 0.0, 0.0, 0.0};
    const GQCP::Nucleus o {8, 1.0, 0.0, 0.0};
    const GQCP::Nucleus h2 {1, 2.0, 0.0, 0.0};
    const GQCP::Molecule molecule {{h1, o, h2}};

    const auto B = GQCP::HomogeneousMagneticField {{0.0, 0.0, 0.0}};  // Gauge origin at cartesian origin.
    // const auto B = GQCP::HomogeneousMagneticField {{0.0, 0.0, 1.0}, {0.4, 12.2, -0.789}};  // Gauge origin at random point in space.
    auto spin_orbital_basis = GQCP::RSpinOrbitalBasis<GQCP::complex, GQCP::LondonGTOShell> {molecule, "STO-3G", B};
    // auto spin_orbital_basis = GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> {molecule, "STO-3G"};
    
    const GQCP::Vector<double, 3> r = {1.3, -0.9, 3.7};

    // get ao values
    const auto ao_values = spin_orbital_basis.evalBasisSetAtPoint(r);

    // print output
    std::cout << "AO values at point r = (" << r.transpose() << "):" << std::endl;
    for (size_t i = 0; i < ao_values.size(); ++i) {
        std::cout << "AO[" << i << "] = " << ao_values[i] << std::endl;
    }

    return 0;
}