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
#include "Utilities/aliases.hpp"
#include "Utilities/literals.hpp"
#include "gqcpy/include/interfaces.hpp"

#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Bind a specific templated second-quantized Hamiltonian.
 * 
 *  @tparam Hamiltonian         The type of the second-quantized Hamiltonian.
 *  
 *  @param name                 The final class name of the second-quantized Hamiltonian in the Python bindings.
 *  @param description          The description of the Python class.
 */
template <typename Hamiltonian>
py::class_<Hamiltonian> bindSQHamiltonian(py::module& module, const std::string& name, const std::string& description) {

    // Alias some types related to the given Hamiltonian.
    using Transformation = typename Hamiltonian::Transformation;
    using ScalarSQOneElectronOperator = typename Hamiltonian::ScalarSQOneElectronOperator;
    using ScalarSQTwoElectronOperator = typename Hamiltonian::ScalarSQTwoElectronOperator;


    // Define the Python class related to the Hamiltonian.
    py::class_<Hamiltonian> py_Hamiltonian {module, name.c_str(), description.c_str()};


    // Expose the actual Python bindings.
    py_Hamiltonian

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
            "coreContributions",
            [](const Hamiltonian& hamiltonian) {
                return hamiltonian.coreContributions();
            },
            "Return a read-only reference to the contributions to the total one-electron interaction operator.")

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


    return py_Hamiltonian;
}


/**
 *  Bind all types of `SQHamiltonian`s.
 */
void bindSQHamiltonians(py::module& module) {

    auto py_RSQHamiltonian_d = bindSQHamiltonian<RSQHamiltonian<double>>(module, "RSQHamiltonian_d", "A (real) second-quantized Hamiltonian expressed in a restricted spin-orbital basis.");
    py_RSQHamiltonian_d
        .def_static(
            "FromHubbard",
            [](const HubbardHamiltonian<double>& hubbard_hamiltonian) {
                return RSQHamiltonian<double>::FromHubbard(hubbard_hamiltonian);
            });


    bindSQHamiltonian<RSQHamiltonian<complex>>(module, "RSQHamiltonian_cd", "A (complex) second-quantized Hamiltonian expressed in a restricted spin-orbital basis.");

    auto py_USQHamiltonian_d = bindSQHamiltonian<USQHamiltonian<double>>(module, "USQHamiltonian_d", "A (real) second-quantized Hamiltonian expressed in an unrestricted spin-orbital basis.");
    py_USQHamiltonian_d
        .def_static(
            "FromRestricted",
            [](const RSQHamiltonian<double>& r_hamiltonian) {
                return USQHamiltonian<double>::FromRestricted(r_hamiltonian);
            });

    bindSQHamiltonian<USQHamiltonian<complex>>(module, "USQHamiltonian_cd", "A (complex) second-quantized Hamiltonian expressed in an unrestricted spin-orbital basis.");

    bindSQHamiltonian<GSQHamiltonian<double>>(module, "GSQHamiltonian_d", "A (real) second-quantized Hamiltonian expressed in a generalized spinor basis.");
    bindSQHamiltonian<GSQHamiltonian<complex>>(module, "GSQHamiltonian_cd", "A (complex) second-quantized Hamiltonian expressed in a generalized spinor basis.");
}


}  // namespace gqcpy
