#include "Basis/ScalarBasis/GTOBasisSet.hpp"
#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Basis/ScalarBasis/ShellSet.hpp"
#include "Mathematical/Functions/CartesianGTO.hpp"
#include "Mathematical/Functions/LondonCartesianGTO.hpp"
#include "Mathematical/Functions/EvaluableLinearCombination.hpp"
#include "Basis/SpinorBasis/RSpinOrbitalBasis.hpp"
#include <array>
#include <iostream>

// -- free helper --
auto evalGradAOsAtPoint(
    const GQCP::RSpinOrbitalBasis<GQCP::complex, GQCP::LondonGTOShell>& spin_basis,
    const GQCP::Vector<double,3>& r
) {
    // get the contracted AOs from the scalar basis
    const auto& AOs = spin_basis.scalarBasis().basisFunctions();
    size_t K = AOs.size();

    // prepare return container: K arrays of 3 components
    std::vector<std::array<GQCP::complex,3>> grad_AO_vals(K);

    for (size_t i = 0; i < K; ++i) {
        const auto& bf     = AOs[i];
        const auto& coeffs = bf.coefficients();   // contraction coefficients
        const auto& prims  = bf.functions();      // primitives

        // accumulate gradient in each direction
        std::array<GQCP::complex,3> sum{0,0,0};
        for (size_t d = 0; d < prims.size(); ++d) {
            auto c = coeffs[d];
            // primitive gradients: a Vector<EvaluableLinearCombination<...>,3>
            auto prim_grad = prims[d].calculatePositionGradient();

            for (int dir = 0; dir < 3; ++dir) {
                // evaluate that linear combination at r
                sum[dir] += c * prim_grad(dir)(r);
            }
        }
        grad_AO_vals[i] = sum;
    }

    return grad_AO_vals;
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
    auto grad_values = evalGradAOsAtPoint(spin_orbital_basis, r);
    std::cout << "AO gradients at r:\n";
    for (size_t i = 0; i < grad_values.size(); ++i) {
        auto& g = grad_values[i];
        std::cout << "  AO["<<i<<"] = (" << g[0] << ", " << g[1] << ", " << g[2] << ")\n";
    }

    return 0;
}