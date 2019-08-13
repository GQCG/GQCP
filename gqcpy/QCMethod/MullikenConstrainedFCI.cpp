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

#include "Molecule/Molecule.hpp"


namespace py = pybind11;


namespace gqcpy {


void bindMolecule(py::module& module) {
    py::class_<GQCP::QCMethod::MullikenConstrainedFCI>(module, "MullikenConstrainedFCI", "A class that represents a collection of nuclei with a number of electrons")
        .def(py::init<const std::vector<GQCP::Nucleus>& , const int>(), py::arg("nuclei"), py::arg("charge") = 0)
        .def("__repr__", [](const GQCP::Molecule& m){ std::ostringstream ss;  ss<<m; return ss.str(); })
        .def_static("ReadXYZ", &GQCP::Molecule::ReadXYZ, "Construct a molecule based on the content of a given .xyz-file.");
}



}  // namespace gqcpy
