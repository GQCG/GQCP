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

#include "ONVBasis/SpinUnresolvedSelectedONVBasis.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Register `SpinUnresolvedSelectedONVBasis` to the gqcpy module and expose parts of its C++ interfaces to Python.
 * 
 *  @param module           The Pybind11 module in which the class should be registered.
 */
void bindSpinUnresolvedSelectedONVBasis(py::module& module) {

    // Define the Python class for `SpinResolvedONVB`.
    py::class_<SpinUnresolvedSelectedONVBasis> py_SpinUnresolvedSelectedONVBasis {module, "SpinUnresolvedSelectedONVBasis", "A spin-resolved ONV basis with a flexible number of (spin-unresolved) ONVs."};

    py_SpinUnresolvedSelectedONVBasis

        /*
         *  MARK: Constructors
         */

        .def(py::init<const size_t, const size_t>(),
             py::arg("M"),
             py::arg("N"),
             "Construct an empty spin-resolved selected ONV basis.")

        .def(py::init<const SpinUnresolvedONVBasis&>(),
             py::arg("onv_basis"),
             "Generate a `SpinUnresolvedSelectedONVBasis` from a spin-unresolved ONV basis.")


        /*
         *  MARK: General information
         */

        .def("numberOfOrbitals",
             &SpinUnresolvedSelectedONVBasis::numberOfOrbitals,
             "Return the number of spinors.")

        .def("numberOfElectrons",
             &SpinUnresolvedSelectedONVBasis::numberOfElectrons,
             "Return the number of electrons, i.e. the number of occupied spinors.")

        .def("dimension",
             &SpinUnresolvedSelectedONVBasis::dimension,
             "Return the dimension of the Fock subspace that is spanned by this selected ONV basis.")


        /*
         *  MARK: Modifying
         */

        .def("expandWith",
             [](SpinUnresolvedSelectedONVBasis& onv_basis, const SpinUnresolvedONV& onv) {
                 onv_basis.expandWith(onv);
             })

        .def("expandWith",
             [](SpinUnresolvedSelectedONVBasis& onv_basis, const std::vector<SpinUnresolvedONV>& onvs) {
                 onv_basis.expandWith(onvs);
             });
}


}  // namespace gqcpy
