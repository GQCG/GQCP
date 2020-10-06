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

#include "Operator/SecondQuantized/SQOneElectronOperator.hpp"
#include "Utilities/aliases.hpp"

#include <pybind11/eigen.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace py = pybind11;

namespace gqcpy {


/**
 *  Since SQOneElectronOperator has a template argument for the representation of the scalar type, we'll have to bind each of them separately. In order to avoid duplicate code, we use a templated binding approach.
 */

/**
 *  Bind a scalar representation of a SQOneElectronOperator.
 * 
 *  @tparam Scalar              the scalar type of the SQOperator
 * 
 *  @param module               the Pybind11 module
 *  @param suffix               the suffix for the gqcpy class name, i.e. "SQOperator" + suffix
 */
template <typename Scalar>
void bindSQOneElectronOperator(py::module& module, const std::string& suffix) {

    // Create a Python binding for ScalarSQOneElectronOperator
    py::class_<GQCP::SQOneElectronOperator<Scalar, 1>>(module,
                                                       ("ScalarSQOneElectronOperator_" + suffix).c_str(),
                                                       "A class that represents a second-quantized one-electron operator")

        // PUBLIC METHODS

        .def(double() * py::self)

        .def(
            "calculateExpectationValue",
            [](const GQCP::SQOneElectronOperator<Scalar, 1>& op, const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& D) {  // use an itermediary Eigen matrix for the Python binding, since Pybind11 doesn't accept our types that are derived from Eigen::Matrix
                return op.calculateExpectationValue(GQCP::OneDM<Scalar> {D});
            },
            "Return the expectation value of the scalar one-electron operator given a 1-DM.")

        .def(
            "parameters",
            [](const GQCP::SQOneElectronOperator<Scalar, 1>& op) {
                return op.parameters().Eigen();
            },
            "Return the integrals encapsulated by the second-quantized one-electron operator.")

        .def(
            "rotate",
            [](GQCP::SQOneElectronOperator<Scalar, 1>& sq_one_op, const Eigen::MatrixXd& U) {
                sq_one_op.rotate(GQCP::TransformationMatrix<double> {U});
            },
            "In-place rotate the operator to another basis.",
            py::arg("U"))

        .def(
            "transform",
            [](GQCP::SQOneElectronOperator<Scalar, 1>& sq_one_op, const Eigen::MatrixXd& T) {
                sq_one_op.transform(GQCP::TransformationMatrix<double> {T});
            },
            "In-place transform the operator to another basis.",
            py::arg("T"));


    // Create a Python binding for VectorSQOneElectronOperator
    py::class_<GQCP::SQOneElectronOperator<Scalar, 3>>(module,
                                                       ("VectorSQOneElectronOperator_" + suffix).c_str(),
                                                       "A class that represents a second-quantized one-electron operator with three components.")

        // PUBLIC METHODS

        .def(
            "allParameters",
            [](const GQCP::SQOneElectronOperator<Scalar, 3>& op) {
                const auto all_parameters = op.allParameters();  // returns a std::array<SquareMatrix<Scalar>>

                return all_parameters;
            },
            "Return the integrals encapsulated by the second-quantized one-electron operator.")

        .def(
            "calculateExpectationValue",
            [](const GQCP::SQOneElectronOperator<Scalar, 3>& op, const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& D) {  // use an itermediary Eigen matrix for the Python binding, since Pybind11 doesn't accept our types that are derived from Eigen::Matrix
                return op.calculateExpectationValue(GQCP::OneDM<Scalar> {D});
            },
            "Return the expectation value of the vector one-electron operator given a 1-DM.")

        .def(
            "rotate",
            [](GQCP::SQOneElectronOperator<Scalar, 3>& sq_one_op, const Eigen::MatrixXd& U) {
                sq_one_op.rotate(GQCP::TransformationMatrix<double> {U});
            },
            "In-place rotate the operator to another basis.",
            py::arg("U"))

        .def(
            "transform",
            [](GQCP::SQOneElectronOperator<Scalar, 3>& sq_one_op, const Eigen::MatrixXd& T) {
                sq_one_op.transform(GQCP::TransformationMatrix<double> {T});
            },
            "In-place transform the operator to another basis.",
            py::arg("T"));
}


void bindSQOneElectronOperators(py::module& module) {

    bindSQOneElectronOperator<double>(module, "d");          // suffix 'd' for the class name
    bindSQOneElectronOperator<GQCP::complex>(module, "cd");  // suffix 'cd' for the class name: 'complex double'
}


}  // namespace gqcpy