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


namespace py = pybind11;


namespace gqcpy {


void bindOrbitalSpace(py::module& module) {
    py::class_<GQCP::OrbitalSpace>(module, "OrbitalSpace", "A class that encapsulates occupied, active and virtual orbital indices.")

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
            "Occupied",
            [](const size_t N) {
                return GQCP::OrbitalSpace::Occupied(N);
            },
            py::arg("N"),
            "Create an orbital space with only occupied indices [0, N[.")

        .def_static(
            "OccupiedVirtual",
            [](const size_t N, const size_t M) {
                return GQCP::OrbitalSpace::OccupiedVirtual(N, M);
            },
            py::arg("N"),
            py::arg("M"),
            "Create an orbital space that is separated between occupied and virtual orbitals.")


        // PUBLIC METHODS
        .def(
            "allIndices",
            &GQCP::OrbitalSpace::allIndices,
            "Return all the indices of the spinors.")

        .def(
            "activeIndices",
            &GQCP::OrbitalSpace::activeIndices,
            "Return the active orbital space, i.e. those spinor indices that are occupied in some set of configurations but not in others.")

        .def(
            "description",
            &GQCP::OrbitalSpace::description,
            "Return a textual description of this orbital space.")

        .def(
            "isIndexActive",
            &GQCP::OrbitalSpace::isIndexActive,
            py::arg("p"),
            "Return if the orbital at the given index is in the active orbital space.")

        .def(
            "isIndexOccupied",
            &GQCP::OrbitalSpace::isIndexOccupied,
            py::arg("p"),
            "Return if the orbital at the given index is in the occupied orbital space.")

        .def(
            "isIndexVirtual",
            &GQCP::OrbitalSpace::isIndexVirtual,
            py::arg("p"),
            "Return if the orbital at the given index is in the virtual orbital space.")

        .def(
            "numberOfOrbitals",
            &GQCP::OrbitalSpace::numberOfOrbitals,
            "Return the total number of orbitals (i.e. spatial orbitals or spinors, depending on the context) in this orbital space.")

        .def(
            "numberOfActiveOrbitals",
            &GQCP::OrbitalSpace::numberOfActiveOrbitals,
            "Return the number of active orbitals (i.e. spatial orbitals or spinors, depending on the context) in this orbital space.")

        .def(
            "numberOfOccupiedOrbitals",
            &GQCP::OrbitalSpace::numberOfOccupiedOrbitals,
            "Return the number of occupied orbitals (i.e. spatial orbitals or spinors, depending on the context) in this orbital space.")

        .def(
            "numberOfVirtualOrbitals",
            &GQCP::OrbitalSpace::numberOfVirtualOrbitals,
            "Return the number of virtual orbitals (i.e. spatial orbitals or spinors, depending on the context) in this orbital space.")

        .def(
            "occupiedIndices",
            &GQCP::OrbitalSpace::occupiedIndices,
            "Return the occupied orbital space, i.e. the indices of the orbitals that are considered occupied by the electrons.")

        .def(
            "virtualIndices",
            &GQCP::OrbitalSpace::virtualIndices,
            "Return the virtual orbital space, i.e. the indices of the orbitals that are considered virtual by the electrons.");
}


}  // namespace gqcpy
