// This file is part of GQCG-GQCP.
//
// Copyright (C) 2017-2020  the GQCG developers
//
// GQCG-GQCP is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GQCG-GQCP is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-GQCP.  If not, see <http://www.gnu.org/licenses/>.

#include "Basis/Integrals/FunctionalPrimitiveEngine.hpp"
#include "Basis/Integrals/IntegralCalculator.hpp"
#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Basis/ScalarBasis/ShellSet.hpp"
#include "Utilities/aliases.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Register `IntegralCalculator.calculate` to the gqcpy module.
 * 
 *  @param module           The Pybind11 module in which `IntegralCalculator.calculate` should be registered.
 */
void bindIntegralEngine(py::module& module) {

    auto py_module_IntegralCalculator = module.def_submodule("IntegralCalculator");

    py_module_IntegralCalculator.def(
        "calculate",
        [](OneElectronIntegralEngine<FunctionalPrimitiveEngine<double>>& engine, const ShellSet<GTOShell>& left_shell_set, const ShellSet<GTOShell>& right_shell_set) {
            return IntegralCalculator::calculate(engine, left_shell_set, right_shell_set)[0];
        },
        "Calculate all one-electron integrals over the basis functions inside the given shell sets.");
}


}  // namespace gqcpy
