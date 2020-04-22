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
#include "Molecule/Molecule.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace py = pybind11;


namespace gqcpy {


void bindMolecule(py::module& module) {
    py::class_<GQCP::Molecule>(module, "Molecule", "A class that represents a collection of nuclei with a number of electrons")

        .def(py::init<const std::vector<GQCP::Nucleus>&, const int>(),
             py::arg("nuclei"),
             py::arg("charge") = 0)

        .def("__repr__",
             [](const GQCP::Molecule& m) {
                 std::ostringstream ss;
                 ss << m;
                 return ss.str();
             })

        .def_static("HCHain",
                    &GQCP::Molecule::HChain,
                    "Return a H-chain with equal internuclear spacing.")

        .def_static("H2CHain",
                    &GQCP::Molecule::H2Chain,
                    "Return an H2-chain.")

        .def_static("HRingFromDistance",
                    &GQCP::Molecule::HRingFromDistance,
                    "Return a regular H-ring where neighbouring hydrogens are separated by the given distance.")

        .def_static("HRingFromRadius",
                    &GQCP::Molecule::HRingFromRadius,
                    "Return a regular H-ring whose hydrogens are on the circle with the given radius.")

        .def_static("ReadXYZ",
                    &GQCP::Molecule::ReadXYZ,
                    "Construct a molecule based on the content of a given .xyz-file.")


        .def("numberOfElectrons",
             &GQCP::Molecule::numberOfElectrons,
             "Return the number of electrons in the molecule.")

        .def("numberOfElectronPairs",
             &GQCP::Molecule::numberOfElectronPairs,
             "Return the number of electron pairs in this molecule. For odd numbers of electrons, the number of electron pairs is equal to that of the (N-1)-even-electron system.");
}


}  // namespace gqcpy
