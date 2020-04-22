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
#include "QCModel/Geminals/AP1roGGeminalCoefficients.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace py = pybind11;


namespace gqcpy {


void bindAP1roGGeminalCoefficients(py::module& module) {

    py::class_<GQCP::AP1roGGeminalCoefficients>(module,
                                                "AP1roGGeminalCoefficients",
                                                "Geminal coefficients for an AP1roG wave function.")

        // Define an initializer through an Eigen::Matrix.
        .def(py::init(
                 [](const Eigen::MatrixXd& G) {
                     return GQCP::AP1roGGeminalCoefficients(GQCP::MatrixX<double> {G});
                 }),
             py::arg("G"))

        .def(py::init<size_t, size_t>(),
             py::arg("N_P"),
             py::arg("K"))

        .def_static("WeakInteractionLimit",
                    &GQCP::AP1roGGeminalCoefficients::WeakInteractionLimit,
                    py::arg("sq_hamiltonian"),
                    py::arg("N_P"),
                    "Return the AP1roG geminal coefficients in the weak interaction limit.")

        .def("asMatrix",
             &GQCP::AP1roGGeminalCoefficients::asMatrix,
             "Return the total geminal coefficient matrix, including the identity matrix block.");
}


}  // namespace gqcpy
