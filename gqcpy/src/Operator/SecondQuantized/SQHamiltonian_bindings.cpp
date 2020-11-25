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

#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "gqcpy/include/interfaces.hpp"

#include <pybind11/operators.h>
#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Bind a specific templated second-quantized Hamiltonian.
 * 
 *  @tparam Hamiltonian         The type of the second-quantized Hamiltonian.
 *  @tparam SpinorBasis         The type of spinor basis that is related to the Hamiltonian. TODO: This can be rewritten after #423.
 *  
 *  @param name                 The final class name of the second-quantized Hamiltonian in the Python bindings.
 *  @param description          The description of the Python class.
 */
template <typename Hamiltonian, typename SpinorBasis>
void bindSQHamiltonian(py::module& module, const std::string& name, const std::string& description) {

    // Alias some types related to the given Hamiltonian.
    using Transformation = typename Hamiltonian::Transformation;
    using ScalarSQOneElectronOperator = typename Hamiltonian::ScalarSQOneElectronOperator;
    using ScalarSQTwoElectronOperator = typename Hamiltonian::ScalarSQTwoElectronOperator;


    // Define the Python class related to the Hamiltonian.
    py::class_<Hamiltonian> py_Hamiltonian {module, name.c_str(), description.c_str()};


    // Expose the actual Python bindings.
    py_Hamiltonian

        /*
         *  MARK: Named constructors
         */

        .def_static(
            "Molecular",
            [](const SpinorBasis& spinor_basis, const Molecule& molecule) {
                return Hamiltonian::Molecular(spinor_basis, molecule);
            },
            py::arg("spinor_basis"),
            py::arg("molecule"),
            "Construct the molecular Hamiltonian in a spinor basis.")


        /*
         *  MARK: Access
         */

        .def(
            "core",
            [](const Hamiltonian& hamiltonian) {
                return hamiltonian.core();
            },
            "Return a read-only reference to the total one-electron interaction operator, i.e. the core Hamiltonian.")

        .def(
            "twoElectron",
            [](const Hamiltonian& hamiltonian) {
                return hamiltonian.twoElectron();
            },
            "Return a read-only reference to the total one-electron interaction operator, i.e. the core Hamiltonian.")


        /*
         *  MARK: Operations related to one-electron operators
         */

        .def(py::self += ScalarSQOneElectronOperator())

        .def(py::self + ScalarSQOneElectronOperator())

        .def(ScalarSQOneElectronOperator() + py::self)

        .def(py::self -= ScalarSQOneElectronOperator())

        .def(py::self - ScalarSQOneElectronOperator())


        /*
         *  MARK: Operations related to two-electron operators
         */

        .def(py::self += ScalarSQTwoElectronOperator())

        .def(py::self + ScalarSQTwoElectronOperator())

        .def(ScalarSQTwoElectronOperator() + py::self)

        .def(py::self -= ScalarSQTwoElectronOperator())

        .def(py::self - ScalarSQTwoElectronOperator());


    // Expose the `BasisTransformable` interface.
    bindBasisTransformableInterface(py_Hamiltonian);
}


/**
 *  Bind all types of `SQHamiltonian`s.
 */
void bindSQHamiltonians(py::module& module) {

    bindSQHamiltonian<RSQHamiltonian<double>, RSpinOrbitalBasis<double, GTOShell>>(module, "RSQHamiltonian", "A second-quantized Hamiltonian expressed in a restricted spin-orbital basis.");
    bindSQHamiltonian<USQHamiltonian<double>, USpinOrbitalBasis<double, GTOShell>>(module, "USQHamiltonian", "A second-quantized Hamiltonian expressed in an unrestricted spin-orbital basis.");
    bindSQHamiltonian<GSQHamiltonian<double>, GSpinorBasis<double, GTOShell>>(module, "GSQHamiltonian", "A second-quantized Hamiltonian expressed in a generalized spinor basis.");
}


}  // namespace gqcpy
