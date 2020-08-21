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

#include "Operator/SecondQuantized/SQTwoElectronOperator.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace py = pybind11;


namespace gqcpy {


/**
 *  Convert a rank-four Eigen::Tensor as a numpy array
 * 
 *  @param tensor       the tensor that should be converted to the numpy array
 * 
 *  @return the corresponding numpy array
 */
template <typename T>
py::array_t<T> asNumpyArray(const Eigen::Tensor<T, 4>& tensor) {

    // Implementation adapted from https://github.com/pybind/pybind11/issues/1377
    const auto shape = tensor.dimensions();

    return py::array_t<T>(shape,
                          {shape[0] * shape[1] * shape[2] * sizeof(T), shape[1] * shape[2] * sizeof(T), shape[2] * sizeof(T), sizeof(T)},  // strides
                          tensor.data()                                                                                                    // data pointer
    );
}


void bindSQTwoElectronOperator(py::module& module) {
    py::class_<GQCP::SQTwoElectronOperator<double, 1>>(module, "SQTwoElectronOperator", "A class that represents a real, second-quantized two-electron operator")

        // PUBLIC METHODS

        .def(
            "calculateExpectationValue",
            [](const GQCP::SQTwoElectronOperator<double, 1>& op, const Eigen::Tensor<double, 4>& d) {  // use an itermediary Eigen Tensor for the Python binding, since Pybind11 doesn't accept our types that are derived from Eigen::Tensor
                return op.calculateExpectationValue(GQCP::TwoDM<double> {d});
            },
            "Return the expectation value of the scalar two-electron operator given a 2-DM.")

        .def(
            "parameters",
            [](const GQCP::SQTwoElectronOperator<double, 1>& op) {
                return asNumpyArray(op.parameters().Eigen());
            },
            "Return the integrals encapsulated by the second-quantized two-electron operator")

        .def(
            "rotate",
            [](GQCP::SQTwoElectronOperator<double, 1>& sq_two_op, const Eigen::MatrixXd& U) {
                sq_two_op.rotate(GQCP::TransformationMatrix<double> {U});
            },
            "In-place rotate the operator to another basis.",
            py::arg("U"))

        .def(
            "transform",
            [](GQCP::SQTwoElectronOperator<double, 1>& sq_two_op, const Eigen::MatrixXd& T) {
                sq_two_op.transform(GQCP::TransformationMatrix<double> {T});
            },
            "In-place transform the operator to another basis.",
            py::arg("T"));
}


}  // namespace gqcpy
