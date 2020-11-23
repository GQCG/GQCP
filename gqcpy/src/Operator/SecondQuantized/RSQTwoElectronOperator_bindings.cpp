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

#include "Operator/SecondQuantized/RSQTwoElectronOperator.hpp"
#include "gqcpy/include/utilities.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace py = pybind11;


namespace gqcpy {

using namespace GQCP;


void bindRSQTwoElectronOperator(py::module& module) {
    py::class_<ScalarRSQTwoElectronOperator<double>>(module, "ScalarRSQTwoElectronOperator", "A class that represents a real, scalar, second-quantized two-electron operator.")

        // PUBLIC METHODS

        // .def(
        //     "calculateExpectationValue",
        //     [](const ScalarRSQTwoElectronOperator<double>& op, const Eigen::Tensor<double, 4>& d) {  // use an itermediary Eigen Tensor for the Python binding, since Pybind11 doesn't accept our types that are derived from Eigen::Tensor
        //         return op.calculateExpectationValue(TwoDM<double> {d});
        //     },
        //     "Return the expectation value of the scalar two-electron operator given a 2-DM.")

        .def(
            "parameters",
            [](const ScalarRSQTwoElectronOperator<double>& op) {
                return asNumpyArray(op.parameters().Eigen());
            },
            "Return the integrals encapsulated by the second-quantized two-electron operator.");
}


}  // namespace gqcpy
