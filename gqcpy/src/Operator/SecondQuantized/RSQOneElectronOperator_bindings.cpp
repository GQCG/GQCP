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

#include "Operator/SecondQuantized/RSQOneElectronOperator.hpp"
#include "Utilities/aliases.hpp"

#include <pybind11/eigen.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace py = pybind11;

namespace gqcpy {

using namespace GQCP;


/**
 *  Since RSQOneElectronOperator has a template argument for the representation of the scalar type, we'll have to bind each of them separately. In order to avoid duplicate code, we use a templated binding approach.
 */

/**
 *  Bind a scalar representation of an `RSQOneElectronOperator`.
 * 
 *  @tparam Scalar              The scalar type used for a single parameter: real or complex.
 * 
 *  @param module               The Pybind11 module.
 *  @param suffix               The suffix for the gqcpy class name, i.e. "RSQOneElectronOperator_" + suffix.
 */
template <typename Scalar>
void bindRSQOneElectronOperator(py::module& module, const std::string& suffix) {

    // Create a Python binding for ScalarRSQOneElectronOperator.
    py::class_<ScalarRSQOneElectronOperator<Scalar>>(module,
                                                     ("ScalarRSQOneElectronOperator_" + suffix).c_str(),
                                                     "A restricted (scalar) one-electron operator, which is suited for expressing non-relativistic (spin-free) one-electron operators.")

        /*
         *  MARK: Calculations
         */
        // .def(
        //     "calculateExpectationValue",
        //     &ScalarRSQOneElectronOperator<Scalar>::calculateExpectationValue,
        //     "Return the expectation value of all components of the one-electron operator")


        /*
         *  MARK: Conforming to `VectorSpaceArithmetic`
         */
        .def(double() * py::self)


        /*
         *  MARK: Parameter access
         */

        .def(
            "parameters",
            [](const ScalarRSQOneElectronOperator<Scalar>& op) {
                return op.parameters().Eigen();
            },
            "A read-only matrix representation of the parameters/matrix elements/integrals of one of the tensor components of this operator.");


    // Create a Python binding for VectorSQOneElectronOperator.
    py::class_<VectorRSQOneElectronOperator<Scalar>>(module,
                                                     ("VectorRSQOneElectronOperator_" + suffix).c_str(),
                                                     "A restricted one-electron operator with three tensor components.")

        // PUBLIC METHODS

        .def(
            "allParameters",
            [](const VectorRSQOneElectronOperator<Scalar>& op) {
                return op.allParameters();
            },
            "Return a vector of read-only matrix representions of the parameters/matrix elements/integrals of this operator.")

        // .def(
        //     "calculateExpectationValue",
        //     &ScalarRSQOneElectronOperator<Scalar>::calculateExpectationValue,
        //     "Return the expectation value of all components of the one-electron operator")
        ;
}


void bindRSQOneElectronOperators(py::module& module) {

    bindRSQOneElectronOperator<double>(module, "d");   // Suffix 'd' for the class name.
    bindRSQOneElectronOperator<complex>(module, "c");  // Suffix 'cd' for the class name: complex.
}


}  // namespace gqcpy