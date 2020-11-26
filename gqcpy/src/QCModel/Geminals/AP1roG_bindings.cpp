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

#include "QCModel/Geminals/AP1roG.hpp"

#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


void bindQCModelAP1roG(py::module& module) {
    py::class_<QCModel::AP1roG>(module, "QCModel_AP1roG", "The AP1roG wave function model.")

        // PUBLIC METHODS

        .def("geminalCoefficients",
             &QCModel::AP1roG::geminalCoefficients,
             "Return the corresponding geminal coefficients of these AP1roG model parameters.");
}


}  // namespace gqcpy
