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
#include "version.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


void bindVersion(py::module& module) {
    py::class_<GQCP::Version>(module, "Version", "Information on the GQCP version")
        .def(
            "major",
            &Version::major,
            "The GQCP major version.")

        .def(
            "minor",
            &Version::minor,
            "The GQCP minor version.")

        .def(
            "patch",
            &Version::patch,
            "The GQCP patch version.")

        .def(
            "full",
            &Version::full,
            "The full GQCP version.")

        .def(
            "git_sha1",
            &Version::git_SHA1,
            "The current GQCP commit SHA1.");
}


}  // namespace gqcpy
