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

#include "Operator/FirstQuantized/Operator.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


void bindOperator(py::module& module) {

    py::class_<NuclearDipoleOperator>(module, "NuclearDipoleOperator", "The nuclear dipole operator.")

        // PUBLIC METHODS

        .def("value",
             &NuclearDipoleOperator::value,
             "Return the value of this nuclear dipole operator.");


    py::class_<NuclearRepulsionOperator>(module, "NuclearRepulsionOperator", "The nuclear repulsion operator.")

        // PUBLIC METHODS

        .def("value",
             &NuclearRepulsionOperator::value,
             "Return the scalar value of this nuclear repulsion operator.");


    py::class_<Operator>(module, "Operator", "A class that is used to construct operators using static methods, much like a factory class.")

        // PUBLIC METHODS

        .def_static(
            "NuclearDipole",
            [](const Molecule& molecule, const Eigen::Vector3d& o) {
                return Operator::NuclearDipole(molecule, Vector<double, 3>(o));
            },
            py::arg("molecule"),
            py::arg("o") = Eigen::Vector3d::Zero(),
            "Return a NuclearDipoleOperator.")

        .def_static("NuclearRepulsion",
                    &Operator::NuclearRepulsion,
                    "Return a NuclearRepulsionOperator.");
}


}  // namespace gqcpy
