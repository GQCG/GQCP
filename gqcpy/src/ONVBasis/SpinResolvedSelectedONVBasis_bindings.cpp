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

#include "ONVBasis/SpinResolvedSelectedONVBasis.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Register `SpinResolvedSelectedONVBasis` to the gqcpy module and expose parts of its C++ interfaces to Python.
 * 
 *  @param module           The Pybind11 module in which the class should be registered.
 */
void bindSpinResolvedSelectedONVBasis(py::module& module) {

    // Define the Python class for `SpinResolvedONVB`.
    py::class_<SpinResolvedSelectedONVBasis> py_SpinResolvedSelectedONVBasis {module, "SpinResolvedSelectedONVBasis", "A spin-resolved ONV basis with a flexible number of (spin-resolved) ONVs."};

    py_SpinResolvedSelectedONVBasis

        /*
         *  MARK: Constructors
         */

        .def(py::init<const size_t, const size_t, const size_t>(),
             py::arg("K"),
             py::arg("N_alpha"),
             py::arg("N_beta"),
             "Construct an empty spin-resolved selected ONV basis.")

        .def(py::init<const SeniorityZeroONVBasis&>(),
             py::arg("onv_basis"),
             "Generate a `SpinResolvedSelectedONVBasis` from a seniority-zero ONV basis.")

        .def(py::init<const SpinResolvedONVBasis&>(),
             py::arg("onv_basis"),
             "Generate a `SpinResolvedSelectedONVBasis` from a full spin-resolved ONV basis.")


        /*
         *  MARK: Named constructors
         */

        .def_static(
            "CIS",
            &SpinResolvedSelectedONVBasis::CIS,
            "Return a CI singles-equivalent `SpinResolvedSelectedONVBasis`.")


        /*
         *  MARK: General information
         */

        .def("numberOfOrbitals",
             &SpinResolvedSelectedONVBasis::numberOfOrbitals,
             "Return the number of spin-orbitals (equal for alpha and beta).")

        .def("numberOfAlphaElectrons",
             &SpinResolvedSelectedONVBasis::numberOfAlphaElectrons,
             "Return the number of alpha electrons, i.e. the number of occupied alpha spin-orbitals.")

        .def("numberOfBetaElectrons",
             &SpinResolvedSelectedONVBasis::numberOfBetaElectrons,
             "Return the number of beta electrons, i.e. the number of occupied beta spin-orbitals.")

        .def("dimension",
             &SpinResolvedSelectedONVBasis::dimension,
             "Return the dimension of the Fock subspace that is spanned by this selected ONV basis.")

        .def("onvWithIndex",
             &SpinResolvedSelectedONVBasis::onvWithIndex,
             py::arg("i"),
             "The ONV that corresponds to the given index/address.")

        /*
         *  MARK: Modifying
         */

        .def("expandWith",
             [](SpinResolvedSelectedONVBasis& onv_basis, const SpinResolvedONV& onv) {
                 onv_basis.expandWith(onv);
             })

        .def("expandWith",
             [](SpinResolvedSelectedONVBasis& onv_basis, const std::vector<SpinResolvedONV>& onvs) {
                 onv_basis.expandWith(onvs);
             });
}


}  // namespace gqcpy
