// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2019  the GQCG developers
// 
// GQCG-gqcp is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// GQCG-gqcp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-gqcp.  If not, see <http://www.gnu.org/licenses/>.
// 
#include "Operator/SecondQuantized/SQOneElectronOperator.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace py = pybind11;

namespace gqcpy {


void bindSQOneElectronOperator(py::module& module) {

    // Create a Python binding for ScalarSQOneElectronOperator
    py::class_<GQCP::SQOneElectronOperator<double, 1>>(module, "ScalarSQOneElectronOperator", "A class that represents a real, second-quantized one-electron operator")

        .def("calculateExpectationValue",
            [ ] (const GQCP::SQOneElectronOperator<double, 1>& op, const Eigen::MatrixXd& D) {  // use an itermediary Eigen matrix for the Python binding, since Pybind11 doesn't accept our types that are derived from Eigen::Matrix
                return op.calculateExpectationValue(GQCP::OneRDM<double>{D});
            },
            "Return the expectation value of the scalar one-electron operator given a 1-DM."
        )

        .def("parameters", 
            [ ] (const GQCP::SQOneElectronOperator<double, 1>& op) {
                return op.parameters().Eigen();
            },
            "Return the integrals encapsulated by the second-quantized one-electron operator."
        )
    ;


    // Create a Python binding for VectorSQOneElectronOperator
    py::class_<GQCP::VectorSQOneElectronOperator<double>>(module, "VectorSQOneElectronOperator", "A class that represents a real, second-quantized one-electron operator with three components.")

        .def("allParameters",
            [ ] (const GQCP::VectorSQOneElectronOperator<double>& op) {
                const auto all_parameters = op.allParameters();  // returns a std::array<QCMatrix<double>>
                
                return all_parameters;
            },
            "Return the integrals encapsulated by the second-quantized one-electron operator."
        )


        .def("calculateExpectationValue",
            [ ] (const GQCP::VectorSQOneElectronOperator<double>& op, const Eigen::MatrixXd& D) {  // use an itermediary Eigen matrix for the Python binding, since Pybind11 doesn't accept our types that are derived from Eigen::Matrix
                return op.calculateExpectationValue(GQCP::OneRDM<double>{D});
            },
            "Return the expectation value of the vector one-electron operator given a 1-DM."
        )
    ;
}


}  // namespace gqcpy