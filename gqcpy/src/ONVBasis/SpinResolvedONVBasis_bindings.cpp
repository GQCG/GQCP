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

#include "ONVBasis/SpinResolvedONVBasis.hpp"
#include "gqcpy/include/interfaces.hpp"

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


void bindSpinResolvedONVBasis(py::module& module) {
    py::class_<SpinResolvedONVBasis> py_SpinResolvedONVBasis {module, "SpinResolvedONVBasis", "A full spin-resolved ONV basis."};


    py_SpinResolvedONVBasis

        // CONSTRUCTORS

        .def(py::init<const size_t, const size_t, const size_t>(),
             py::arg("K"),
             py::arg("N_alpha"),
             py::arg("N_beta"))


        // PUBLIC METHODS
        .def(
            "compoundAddress",
            &SpinResolvedONVBasis::compoundAddress,
            py::arg("I_alpha"),
            py::arg("I_beta"),
            "Calculate the compound address of an ONV represented by the two given alpha- and beta-addresses.")

        .def("dimension",
             &SpinResolvedONVBasis::dimension)

        .def(
            "forEach",
            &SpinResolvedONVBasis::forEach,
            py::arg("callback"),
            "Iterate over all ONVs, and apply the given callback function.");


    // Expose the `SpinResolvedBase` interface to the Python class.
    bindSpinResolvedBaseInterface(py_SpinResolvedONVBasis);
}


}  // namespace gqcpy
