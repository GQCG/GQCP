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

#include "ONVBasis/SeniorityZeroONVBasis.hpp"
#include "ONVBasis/SpinResolvedFrozenONVBasis.hpp"
#include "ONVBasis/SpinResolvedONVBasis.hpp"
#include "ONVBasis/SpinResolvedSelectedONVBasis.hpp"
#include "QCModel/CI/LinearExpansion.hpp"

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>


namespace py = pybind11;


namespace gqcpy {


/**
 *  Since LinearExpansion has a template argument for the representation of the ONVBasis, we'll have to bind each of them separately. In order to avoid duplicate code, we use a templated binding approach.
 */

/**
 *  Bind a templated LinearExpansion class.
 * 
 *  @tparam ONVBasis            the type of the ONV basis
 * 
 *  @param module               the Pybind11 module
 *  @param suffix               the suffix for the gqcpy class name, i.e. "LinearExpansion" + suffix
 *  @param description          the description for the gqcpy class
 */
template <typename ONVBasis>
void bindLinearExpansion(py::module& module, const std::string& suffix, const std::string& description) {
    py::class_<GQCP::LinearExpansion<ONVBasis>>(module,
                                                ("LinearExpansion_" + suffix).c_str(),
                                                description.c_str())

        // PUBLIC METHODS
        .def("coefficients",
             &GQCP::LinearExpansion<ONVBasis>::coefficients,
             "Return the expansion coefficients of this linear expansion wave function model.");
}


/**
 *  A template specialization for the binding of GQCP::LinearExpansion<GQCP::SeniorityZeroONVBasis>, because it has an additional function to be bound.
 * 
 *  @param module               the Pybind11 module
 *  @param suffix               the suffix for the gqcpy class name, i.e. "LinearExpansion" + suffix
 *  @param description          the description for the gqcpy class
 */
template <>
void bindLinearExpansion<GQCP::SeniorityZeroONVBasis>(py::module& module, const std::string& suffix, const std::string& description) {

    py::class_<GQCP::LinearExpansion<GQCP::SeniorityZeroONVBasis>>(module,
                                                                   ("LinearExpansion_" + suffix).c_str(),
                                                                   description.c_str())

        // PUBLIC METHODS
        .def(
            "calculate1DM",
            [](const GQCP::LinearExpansion<GQCP::SeniorityZeroONVBasis>& linear_expansion) {
                return linear_expansion.calculate1DM();
            },
            "Return the one-electron density matrix (1-DM) for a seniority-zero wave function expansion.")

        .def("coefficients",
             &GQCP::LinearExpansion<GQCP::SeniorityZeroONVBasis>::coefficients,
             "Return the expansion coefficients of this linear expansion wave function model.");
}


/**
 *  A template specialization for the binding of GQCP::LinearExpansion<GQCP::SpinResolvedONVBasis>, because it has an additional function to be bound.
 * 
 *  @param module               the Pybind11 module
 *  @param suffix               the suffix for the gqcpy class name, i.e. "LinearExpansion" + suffix
 *  @param description          the description for the gqcpy class
 */
template <>
void bindLinearExpansion<GQCP::SpinResolvedONVBasis>(py::module& module, const std::string& suffix, const std::string& description) {

    py::class_<GQCP::LinearExpansion<GQCP::SpinResolvedONVBasis>>(module,
                                                                  ("LinearExpansion_" + suffix).c_str(),
                                                                  description.c_str())

        // CONSTRUCTORS
        .def_static(
            "FromONVProjection",
            [](const GQCP::SpinResolvedONV& onv, const GQCP::RSpinorBasis<double, GQCP::GTOShell>& r_spinor_basis, const GQCP::USpinorBasis<double, GQCP::GTOShell>& u_spinor_basis) {
                return GQCP::LinearExpansion<GQCP::SpinResolvedONVBasis>::FromONVProjection(onv, r_spinor_basis, u_spinor_basis);
            },
            py::arg("onv"),
            py::arg("r_spinor_basis"),
            py::arg("u_spinor_basis"),
            "Create the linear expansion of the given spin-resolved ONV that is expressed in the given USpinorBasis, by projection onto the spin-resolved ONVs expressed with respect to the given RSpinorBasis.")

        // PUBLIC METHODS
        .def("coefficients",
             &GQCP::LinearExpansion<GQCP::SpinResolvedONVBasis>::coefficients,
             "Return the expansion coefficients of this linear expansion wave function model.")

        .def(
            "forEach",
            [](const GQCP::LinearExpansion<GQCP::SpinResolvedONVBasis>& linear_expansion, const std::function<void(const double, const GQCP::SpinResolvedONV)>& callback) {
                return linear_expansion.forEach(callback);
            },
            py::arg("callback"),
            "Iterate over all expansion coefficients and corresponding ONVs, and apply the given callback function.");
}


/**
 *  A template specialization for the binding of GQCP::LinearExpansion<GQCP::SpinUnresolvedONVBasis>, because it has an additional function to be bound.
 * 
 *  @param module               the Pybind11 module
 *  @param suffix               the suffix for the gqcpy class name, i.e. "LinearExpansion" + suffix
 *  @param description          the description for the gqcpy class
 */
template <>
void bindLinearExpansion<GQCP::SpinUnresolvedONVBasis>(py::module& module, const std::string& suffix, const std::string& description) {

    py::class_<GQCP::LinearExpansion<GQCP::SpinUnresolvedONVBasis>>(module,
                                                                    ("LinearExpansion_" + suffix).c_str(),
                                                                    description.c_str())

        // CONSTRUCTORS
        .def_static(
            "FromONVProjection",
            [](const GQCP::SpinUnresolvedONV& onv_of, const GQCP::GSpinorBasis<double, GQCP::GTOShell>& spinor_basis_on, const GQCP::GSpinorBasis<double, GQCP::GTOShell>& spinor_basis_of) {
                return GQCP::LinearExpansion<GQCP::SpinUnresolvedONVBasis>::FromONVProjection(onv_of, spinor_basis_on, spinor_basis_of);
            },
            py::arg("onv_of"),
            py::arg("spinor_basis_on"),
            py::arg("spinor_basis_of"),
            "Create the linear expansion of the given spin-unresolved ONV that is expressed in the given GSpinorBasis, by projection onto the spin-resolved ONVs expressed with respect to another given GSpinorBasis.")

        // PUBLIC METHODS
        .def("coefficients",
             &GQCP::LinearExpansion<GQCP::SpinUnresolvedONVBasis>::coefficients,
             "Return the expansion coefficients of this linear expansion wave function model.");
}


void bindLinearExpansions(py::module& module) {

    bindLinearExpansion<GQCP::SeniorityZeroONVBasis>(module, "SeniorityZero", "The linear expansion (configuration interaction) wave function model in a seniority-zero ONV basis.");
    bindLinearExpansion<GQCP::SpinResolvedFrozenONVBasis>(module, "SpinResolvedFrozen", "The linear expansion (configuration interaction) wave function model in a frozen core spin-resolved ONV basis.");
    bindLinearExpansion<GQCP::SpinResolvedONVBasis>(module, "SpinResolved", "The linear expansion (configuration interaction) wave function model in a spin-resolved ONV basis.");
    bindLinearExpansion<GQCP::SpinResolvedSelectedONVBasis>(module, "SpinResolvedSelected", "The linear expansion (configuration interaction) wave function model in a spin-resolved selected ONV basis.");
    bindLinearExpansion<GQCP::SpinUnresolvedONVBasis>(module, "SpinUnresolved", "The linear expansion (configuration interaction) wave function model in a spin-unresolved ONV basis.");
}


}  // namespace gqcpy
