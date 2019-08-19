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
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "Molecule/Nucleus.hpp"


namespace py = pybind11;


namespace gqcpy {


void bindNucleus(py::module& module) {
    py::class_<GQCP::Nucleus>(module, "Nucleus", "A class that represents a nucleus: it has a charge and a position in space")
        .def(py::init<const size_t, const double, const double, const double>(), py::arg("Z"), py::arg("x"), py::arg("y"), py::arg("z"));
}



}  // namespace gqcpy
