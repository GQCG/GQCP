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

#include "Operator/FirstQuantized/AngularMomentumOperator.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Register `AngularMomentumOperator` to the gqcpy module and expose parts of its C++ interface to Python.
 * 
 *  @param module           The Pybind11 module in which the class should be registered.
 */
void bindAngularMomentumOperator(py::module& module) {

    py::class_<AngularMomentumOperator> py_AngularMomentumOperator {module, "AngularMomentumOperator", "The (one-electron) angular momentum operator."};

    py_AngularMomentumOperator

        /*
         *  MARK: Constructors
         */

        .def(
            py::init<>(
                [](const Eigen::Vector3d& o) {
                    return AngularMomentumOperator(o);
                }),
            py::arg("o") = Eigen::Vector3d::Zero());
}


}  // namespace gqcpy
