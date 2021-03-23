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

#include "Physical/HomogeneousMagneticField.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Register `HomogeneousMagneticField` to the gqcpy module and expose parts of its C++ interface to Python.
 * 
 *  @param module           The Pybind11 module in which the class should be registered.
 */
void bindHomogeneousMagneticField(py::module& module) {

    py::class_<HomogeneousMagneticField> py_HomogeneousMagneticField {module, "HomogeneousMagneticField", "The first-quantized, molecular electronic Hamiltonian."};

    py_HomogeneousMagneticField

        /*
         *  MARK: Constructors
         */

        .def(py::init<>([](const Eigen::Vector3d& B, const Eigen::Vector3d& G) {
                 return HomogeneousMagneticField(B, G);
             }),
             py::arg("B"),
             py::arg("G") = Eigen::Vector3d::Zero())


        /*
         *  MARK: Physics
         */

        .def("strength",
             &HomogeneousMagneticField::strength,
             "Return the field strength of this homogeneous magnetic field.")

        .def("gaugeOrigin",
             &HomogeneousMagneticField::gaugeOrigin,
             "Return the gauge origin of this homogeneous magnetic field.")

        .def(
            "vectorPotentialAt",
            [](const HomogeneousMagneticField& B, const Eigen::Vector3d& r) {
                return B.vectorPotentialAt(r);
            },
            py::arg("r"),
            "Return the vector potential evaluated at the given point.");
}


}  // namespace gqcpy
