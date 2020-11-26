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

#include "Molecule/Molecule.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


void bindMolecule(py::module& module) {
    py::class_<Molecule>(module, "Molecule", "A class that represents a collection of nuclei with a number of electrons")

        // CONSTRUCTORS
        .def(py::init<const std::vector<Nucleus>&, const int>(),
             py::arg("nuclei"),
             py::arg("charge") = 0)

        .def_static(
            "HChain",
            [](const size_t n, const double spacing, const int charge) {
                return Molecule::HChain(n, spacing, charge);
            },
            py::arg("n"),
            py::arg("spacing"),
            py::arg("charge") = 0,
            "Return a H-chain with equal internuclear spacing.")

        .def_static(
            "H2Chain",
            [](const size_t n, const double a, const double b, const int charge) {
                return Molecule::H2Chain(n, a, b, charge);
            },
            py::arg("n"),
            py::arg("a"),
            py::arg("b"),
            py::arg("charge") = 0,
            "Return an H2-chain.")

        .def_static(
            "HRingFromDistance",
            [](const size_t n, const double distance, const int charge) {
                return Molecule::HRingFromDistance(n, distance, charge);
            },
            py::arg("n"),
            py::arg("distance"),
            py::arg("charge") = 0,
            "Return a regular H-ring where neighbouring hydrogens are separated by the given distance.")

        .def_static(
            "HRingFromRadius",
            [](const size_t n, const double radius, const int charge) {
                return Molecule::HRingFromRadius(n, radius, charge);
            },
            py::arg("n"),
            py::arg("radius"),
            py::arg("charge") = 0,
            "Return a regular H-ring whose hydrogens are on the circle with the given radius.")

        .def_static(
            "ReadXYZ",
            [](const std::string& xyz_filename, const int charge) {
                return Molecule::ReadXYZ(xyz_filename, charge);
            },
            py::arg("xyz_filename"),
            py::arg("charge") = 0,
            "Construct a molecule based on the content of a given .xyz-file.")


        // PUBLIC METHODS

        .def("__str__",
             [](const Molecule& molecule) {
                 return molecule.description();
             })

        .def("numberOfElectrons",
             &Molecule::numberOfElectrons,
             "Return the number of electrons in the molecule.")

        .def("numberOfElectronPairs",
             &Molecule::numberOfElectronPairs,
             "Return the number of electron pairs in this molecule. For odd numbers of electrons, the number of electron pairs is equal to that of the (N-1)-even-electron system.");
}


}  // namespace gqcpy
