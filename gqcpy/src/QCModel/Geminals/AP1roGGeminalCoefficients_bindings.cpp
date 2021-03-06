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

#include "QCModel/Geminals/AP1roGGeminalCoefficients.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


void bindAP1roGGeminalCoefficients(py::module& module) {

    py::class_<AP1roGGeminalCoefficients>(module,
                                          "AP1roGGeminalCoefficients",
                                          "Geminal coefficients for an AP1roG wave function.")

        // CONSTRUCTORS
        .def(py::init(
                 [](const Eigen::MatrixXd& G) {
                     return AP1roGGeminalCoefficients(MatrixX<double> {G});
                 }),
             py::arg("G"))

        .def(py::init<size_t, size_t>(),
             py::arg("N_P"),
             py::arg("K"))

        .def_static("WeakInteractionLimit",
                    &AP1roGGeminalCoefficients::WeakInteractionLimit,
                    py::arg("sq_hamiltonian"),
                    py::arg("N_P"),
                    "Return the AP1roG geminal coefficients in the weak interaction limit.")


        // PUBLIC METHODS

        .def("asMatrix",
             &AP1roGGeminalCoefficients::asMatrix,
             "Return the total geminal coefficient matrix, including the identity matrix block.");
}


}  // namespace gqcpy
