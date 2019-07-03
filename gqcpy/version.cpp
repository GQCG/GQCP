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
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "version.hpp"


namespace py = pybind11;


namespace gqcpy {


void bindVersion(py::module& module) {
    py::class_<GQCP::Version>(module, "Version", "Information on the GQCP version")
        .def_property_readonly_static("major", [] (py::object) { return GQCP::Version::major; }, "The GQCP major version")
        .def_property_readonly_static("minor", [] (py::object) { return GQCP::Version::minor; }, "The GQCP minor version")
        .def_property_readonly_static("patch", [] (py::object) { return GQCP::Version::patch; }, "The GQCP patch version")
        .def_property_readonly_static("full", [] (py::object) { return GQCP::Version::full; }, "The full GQCP version")
        .def_property_readonly_static("git_sha1", [] (py::object) { return GQCP::Version::git_SHA1; }, "The current GQCP commit SHA1")
        ;
}


}  // namespace gqcpy
