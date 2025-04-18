#include "Basis/ScalarBasis/GTOBasisSet.hpp"
#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Basis/ScalarBasis/ShellSet.hpp"
#include "Mathematical/Functions/CartesianGTO.hpp"
#include "Mathematical/Functions/LondonCartesianGTO.hpp"
#include "Mathematical/Functions/EvaluableLinearCombination.hpp"
#include "Basis/SpinorBasis/RSpinOrbitalBasis.hpp"
#include <iostream>

// for determining type of object
#include <typeinfo>
#include <cxxabi.h>


std::vector<std::vector<GQCP::complex, 3>> evalGradBasisSetAtPoint(const GQCP::Vector<double, 3>& r) const {
        // gather basis set AOs from the first spatial orbital (which could eg be a spatial MO)
        // which is expanded in the scalar basis set
        const auto AOs = this->spatialOrbitals()[0].functions();
        // init vector in which to gather each AO's value at r
        std::vector<std::vector<GQCP::complex, 3>> grad_AO_vals;
        grad_AO_vals.reserve(this->numberOfSpatialOrbitals()); //n_AO = n_MO
        std::vector<Vector<EvaluableLinearCombination<GQCP::complex, EvaluableLinearCombination<GQCP::complex, LondonCartesianGTO>>, 3>> basis_function_gradients {this->numberOfSpatialOrbitals()};

        // loop through AOs
        for (size_t i= 0; i < this->numberOfSpatialOrbitals(); i++) {
            const auto& basis_function = AOs[i];
            // loop through its primitives
            const auto contraction_length = basis_function.length();
            const auto& contraction_coefficients = basis_function.coefficients();
            const auto& primitives = basis_function.functions();

            for (size_t d = 0; d < contraction_length; d++) {
                const auto& contraction_coefficient = contraction_coefficients[d];
                const auto primitive_gradient = primitives[d].calculatePositionGradient();
                grad_AO_vals[i].append(contraction_coefficient, primitive_gradient(m));
            }

            // evaluate value at r, put it in vector
            grad_AO_vals.push_back(grad_AO_vals[i](r));
        }
        
        return AO_vals;
}

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

    // get gradient of LAOs.



    return 0;
}