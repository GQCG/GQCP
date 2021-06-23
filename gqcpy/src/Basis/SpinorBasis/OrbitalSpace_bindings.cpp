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

#include "Basis/SpinorBasis/OrbitalSpace.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Register `OrbitalSpace` to the gqcpy module and expose a part of its C++ interface to Python.
 * 
 *  @param module           The Pybind11 module in which `OrbitalSpace` should be registered.
 */
void bindOrbitalSpace(py::module& module) {
    py::class_<OrbitalSpace>(module, "OrbitalSpace", "A class that encapsulates occupied, active and virtual orbital indices.")

        // CONSTRUCTORS

        .def(py::init<const std::vector<size_t>&, const std::vector<size_t>&, const std::vector<size_t>&>(),
             py::arg("occupied_indices"),
             py::arg("active_indices"),
             py::arg("virtual_indices"),
             "Construct an orbital space from occupied, active and virtual indices.")

        .def(py::init<const std::vector<size_t>&, const std::vector<size_t>&>(),
             py::arg("occupied_indices"),
             py::arg("virtual_indices"),
             "Construct an orbital space only from occupied and virtual indices.")

        .def_static(
            "Implicit",
            [](const std::map<OccupationType, size_t>& counts) {
                return OrbitalSpace::Implicit(counts);
            },
            py::arg("counts"),
            "Create an implicit orbital space with the given dimensions.")


        // PUBLIC METHODS

        .def(
            "description",
            &OrbitalSpace::description,
            "Return a textual description of this orbital space.")

        .def(
            "indices",
            [](const OrbitalSpace& orbital_space) {
                return orbital_space.indices();
            },
            "Return all the indices of the spinors")

        .def(
            "indices",
            [](const OrbitalSpace& orbital_space, const OccupationType type) {
                return orbital_space.indices(type);
            },
            py::arg("type"),
            "Return the indices that belong to the given occupation type.")

        .def(
            "isIndex",
            [](const OrbitalSpace& orbital_space, const OccupationType type, const size_t p) {
                return orbital_space.isIndex(type, p);
            },
            py::arg("type"),
            py::arg("p"),
            "Return if the orbital at the given index is in the given orbital space.")

        .def(
            "numberOfExcitations",
            &OrbitalSpace::numberOfExcitations,
            py::arg("from"),
            py::arg("to"),
            "Return the number of possible excitations between the two orbital sets")

        .def(
            "numberOfOrbitals",
            [](const OrbitalSpace& orbital_space) {
                return orbital_space.numberOfOrbitals();
            },
            "Return the total number of orbitals (i.e. spatial orbitals or spinors, depending on the context) in this orbital space.")

        .def(
            "numberOfOrbitals",
            [](const OrbitalSpace& orbital_space, const OccupationType type) {
                return orbital_space.numberOfOrbitals(type);
            },
            "Return the number of orbitals (i.e. spatial orbitals or spinors, depending on the context) that belong to the given occupation type.");
}


}  // namespace gqcpy
