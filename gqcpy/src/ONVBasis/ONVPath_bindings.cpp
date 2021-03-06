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

#include "ONVBasis/ONVPath.hpp"
#include "ONVBasis/SpinResolvedONVBasis.hpp"
#include "ONVBasis/SpinUnresolvedONVBasis.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


template <typename ONVBasis>
void bindONVPath(py::module& module, const std::string& suffix, const std::string& description) {
    py::class_<ONVPath<ONVBasis>>(module,
                                  ("ONVPath_" + suffix).c_str(),
                                  description.c_str())

        // CONSTRUCTOR

        .def(py::init<const ONVBasis&, const typename ONVBasis::ONV&>(),
             py::arg("onv_basis"),
             py::arg("onv"))


        // PUBLIC METHODS

        .def(
            "address",
            [](const ONVPath<ONVBasis>& path) {
                return path.address();
            },
            "Return the address of the current path.")

        .def(
            "annihilate",
            [](ONVPath<ONVBasis>& path) {
                path.annihilate();
            },
            "According to this path's current state, annihilate the next diagonal arc.")

        .def(
            "annihilate",
            [](ONVPath<ONVBasis>& path, const size_t q, const size_t n) {
                path.annihilate(q, n);
            },
            py::arg("q"),
            py::arg("n"),
            "Annihilate the diagonal arc that starts at coordinate (q, n).")

        .def(
            "create",
            [](ONVPath<ONVBasis>& path) {
                path.create();
            },
            "According to this path's current state, create the next diagonal arc.")

        .def(
            "create",
            [](ONVPath<ONVBasis>& path, const size_t p, const size_t n) {
                path.create(p, n);
            },
            py::arg("p"),
            py::arg("n"),
            "Create the diagonal arc that starts at coordinate (p, n). ")

        .def(
            "electronIndex",
            [](const ONVPath<ONVBasis>& path) {
                return path.orbitalIndex();
            },
            "Return the electron index 'n' that, together with the orbital index 'p' signifies the vertex (p,n) up until which the ONV path construction is finished.")

        .def(
            "isFinished",
            [](const ONVPath<ONVBasis>& path) {
                return path.isFinished();
            },
            "Return If the path's construction is considered finished.")

        .def(
            "isOrbitalIndexValid",
            [](const ONVPath<ONVBasis>& path) {
                return path.isOrbitalIndexValid();
            },
            "Return If the orbital index 'p' is allowed. If p exceeds a certain value, there are not enough electrons for the given set of orbitals.")

        .def(
            "leftTranslateDiagonalArc",
            [](ONVPath<ONVBasis>& path, const size_t p, const size_t n) {
                path.leftTranslateDiagonalArc(p, n);
            },
            py::arg("p"),
            py::arg("q"),
            "Translate the diagonal arc that starts at the coordinate (p,n) to the left, indicating that the current path is 'open' at the vertex (p,n-1) and that the orbital 'p' should be occupied in subsequent path manipulations.")

        .def(
            "leftTranslateDiagonalArcUntilVerticalArc",
            [](ONVPath<ONVBasis>& path) {
                path.leftTranslateDiagonalArcUntilVerticalArc();
            },
            "According to this path's current state, translate diagonal arcs to the left until an unoccupied orbital (vertical arc) is found.")

        .def(
            "orbitalIndex",
            [](const ONVPath<ONVBasis>& path) {
                return path.orbitalIndex();
            },
            "Return the orbital index 'p' that, together with the electron index 'n' signifies the vertex (p,n) up until which the ONV path construction is finished.")

        .def(
            "sign",
            [](const ONVPath<ONVBasis>& path) {
                return path.sign();
            },
            "Return the total phase factor/sign associated to the original path's modification");
}


void bindONVPaths(py::module& module) {

    // bindONVPath<SpinResolvedONVBasis>(module, "SpinResolved", "A path like representation of a spin resolved ONV.");
    bindONVPath<SpinUnresolvedONVBasis>(module, "SpinUnresolved", "A path like representation of a spin unresolved ONV.");
}

}  // namespace gqcpy
