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
#include "ONVBasis/SpinResolvedONVBasis.hpp"
#include "ONVBasis/SpinResolvedSelectedONVBasis.hpp"
#include "QCModel/CI/LinearExpansion.hpp"
#include "gqcpy/include/utilities.hpp"

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


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
    py::class_<LinearExpansion<ONVBasis>>(module,
                                          ("LinearExpansion_" + suffix).c_str(),
                                          description.c_str())

        // CONSTRUCTORS

        .def_static(
            "Constant",
            [](const ONVBasis& onv_basis) {
                return LinearExpansion<ONVBasis>::Constant(onv_basis);
            },
            py::arg("onv_basis"),
            "Return a linear expansion with a normalized coefficient vector (i.e. all the coefficients are equal).")

        .def_static(
            "HartreeFock",
            [](const ONVBasis& onv_basis) {
                return LinearExpansion<ONVBasis>::HartreeFock(onv_basis);
            },
            py::arg("onv_basis"),
            "Return a linear expansion that represents the Hartree-Fock wave function.")

        .def_static(
            "Random",
            [](const ONVBasis& onv_basis) {
                return LinearExpansion<ONVBasis>::Random(onv_basis);
            },
            py::arg("onv_basis"),
            "Return a linear expansion with a random, normalized coefficient vector, with coefficients uniformly distributed in [-1, +1] before any normalization.")


        // PUBLIC METHODS

        .def("coefficients",
             &LinearExpansion<ONVBasis>::coefficients,
             "Return the expansion coefficients of this linear expansion wave function model.");
}


/**
 *  A template specialization for the binding of LinearExpansion<SeniorityZeroONVBasis>, because it has an additional function to be bound.
 * 
 *  @param module               the Pybind11 module
 *  @param suffix               the suffix for the gqcpy class name, i.e. "LinearExpansion" + suffix
 *  @param description          the description for the gqcpy class
 */
template <>
void bindLinearExpansion<SeniorityZeroONVBasis>(py::module& module, const std::string& suffix, const std::string& description) {

    py::class_<LinearExpansion<SeniorityZeroONVBasis>>(module,
                                                       ("LinearExpansion_" + suffix).c_str(),
                                                       description.c_str())

        // CONSTRUCTORS

        .def_static(
            "Constant",
            [](const SeniorityZeroONVBasis& onv_basis) {
                return LinearExpansion<SeniorityZeroONVBasis>::Constant(onv_basis);
            },
            py::arg("onv_basis"),
            "Return a linear expansion with a normalized coefficient vector (i.e. all the coefficients are equal).")

        .def_static(
            "HartreeFock",
            [](const SeniorityZeroONVBasis& onv_basis) {
                return LinearExpansion<SeniorityZeroONVBasis>::HartreeFock(onv_basis);
            },
            py::arg("onv_basis"),
            "Return a linear expansion that represents the Hartree-Fock wave function.")

        .def_static(
            "Random",
            [](const SeniorityZeroONVBasis& onv_basis) {
                return LinearExpansion<SeniorityZeroONVBasis>::Random(onv_basis);
            },
            py::arg("onv_basis"),
            "Return a linear expansion with a random, normalized coefficient vector, with coefficients uniformly distributed in [-1, +1] before any normalization.")


        // PUBLIC METHODS

        .def(
            "calculate1DM",
            [](const LinearExpansion<SeniorityZeroONVBasis>& linear_expansion) {
                return linear_expansion.calculate1DM();
            },
            "Return the one-electron density matrix (1-DM) for a seniority-zero wave function expansion.")

        .def(
            "calculate2DM",
            [](const LinearExpansion<SeniorityZeroONVBasis>& linear_expansion) {
                return asNumpyArray(linear_expansion.calculate2DM().Eigen());
            },
            "Return the two-electron density matrix (2-DM) for a seniority-zero wave function expansion.")

        .def("coefficients",
             &LinearExpansion<SeniorityZeroONVBasis>::coefficients,
             "Return the expansion coefficients of this linear expansion wave function model.")

        .def(
            "calculateSpinResolved1DM",
            [](const LinearExpansion<SeniorityZeroONVBasis>& linear_expansion) {
                return linear_expansion.calculateSpinResolved1DM();
            },
            "Return the spin-resolved 1-DM.")

        .def(
            "calculateSpinResolved2DM",
            [](const LinearExpansion<SeniorityZeroONVBasis>& linear_expansion) {
                return linear_expansion.calculateSpinResolved2DM();
            },
            "Return the spin resolved 2-DM.");
}


/**
 *  A template specialization for the binding of LinearExpansion<SpinResolvedONVBasis>, because it has an additional function to be bound.
 * 
 *  @param module               the Pybind11 module
 *  @param suffix               the suffix for the gqcpy class name, i.e. "LinearExpansion" + suffix
 *  @param description          the description for the gqcpy class
 */
template <>
void bindLinearExpansion<SpinResolvedONVBasis>(py::module& module, const std::string& suffix, const std::string& description) {

    py::class_<LinearExpansion<SpinResolvedONVBasis>>(module,
                                                      ("LinearExpansion_" + suffix).c_str(),
                                                      description.c_str())

        // CONSTRUCTORS

        .def_static(
            "Constant",
            [](const SpinResolvedONVBasis& onv_basis) {
                return LinearExpansion<SpinResolvedONVBasis>::Constant(onv_basis);
            },
            py::arg("onv_basis"),
            "Return a linear expansion with a normalized coefficient vector (i.e. all the coefficients are equal).")

        .def_static(
            "FromONVProjection",
            [](const SpinResolvedONV& onv, const RSpinOrbitalBasis<double, GTOShell>& r_spinor_basis, const USpinOrbitalBasis<double, GTOShell>& u_spinor_basis) {
                return LinearExpansion<SpinResolvedONVBasis>::FromONVProjection(onv, r_spinor_basis, u_spinor_basis);
            },
            py::arg("onv"),
            py::arg("r_spinor_basis"),
            py::arg("u_spinor_basis"),
            "Create the linear expansion of the given spin-resolved ONV that is expressed in the given USpinOrbitalBasis, by projection onto the spin-resolved ONVs expressed with respect to the given RSpinOrbitalBasis.")

        .def_static(
            "HartreeFock",
            [](const SpinResolvedONVBasis& onv_basis) {
                return LinearExpansion<SpinResolvedONVBasis>::HartreeFock(onv_basis);
            },
            py::arg("onv_basis"),
            "Return a linear expansion that represents the Hartree-Fock wave function.")

        .def_static(
            "Random",
            [](const SpinResolvedONVBasis& onv_basis) {
                return LinearExpansion<SpinResolvedONVBasis>::Random(onv_basis);
            },
            py::arg("onv_basis"),
            "Return a linear expansion with a random, normalized coefficient vector, with coefficients uniformly distributed in [-1, +1] before any normalization.")


        // PUBLIC METHODS

        .def(
            "calculate1DM",
            [](const LinearExpansion<SpinResolvedONVBasis>& linear_expansion) {
                return linear_expansion.calculate1DM();
            },
            "Return the one-electron density matrix (1-DM) for a full spin-resolved wave function expansion.")

        .def(
            "calculate2DM",
            [](const LinearExpansion<SpinResolvedONVBasis>& linear_expansion) {
                return asNumpyArray(linear_expansion.calculate2DM().Eigen());
            },
            "Return the two-electron density matrix (2-DM) for a full spin-resolved wave function expansion.")

        .def("coefficients",
             &LinearExpansion<SpinResolvedONVBasis>::coefficients,
             "Return the expansion coefficients of this linear expansion wave function model.")

        .def(
            "forEach",
            [](const LinearExpansion<SpinResolvedONVBasis>& linear_expansion, const std::function<void(const double, const SpinResolvedONV)>& callback) {
                return linear_expansion.forEach(callback);
            },
            py::arg("callback"),
            "Iterate over all expansion coefficients and corresponding ONVs, and apply the given callback function.")

        .def(
            "calculateSpinResolved1DM",
            [](const LinearExpansion<SpinResolvedONVBasis>& linear_expansion) {
                return linear_expansion.calculateSpinResolved1DM();
            },
            "Return the spin-resolved 1-DM.")

        .def(
            "calculateSpinResolved2DM",
            [](const LinearExpansion<SpinResolvedONVBasis>& linear_expansion) {
                return linear_expansion.calculateSpinResolved2DM();
            },
            "Return the spin resolved 2-DM.");
}


/**
 *  A template specialization for the binding of LinearExpansion<SpinResolvedSelectedONVBasis>, because it has an additional function to be bound.
 * 
 *  @param module               the Pybind11 module
 *  @param suffix               the suffix for the gqcpy class name, i.e. "LinearExpansion" + suffix
 *  @param description          the description for the gqcpy class
 */
template <>
void bindLinearExpansion<SpinResolvedSelectedONVBasis>(py::module& module, const std::string& suffix, const std::string& description) {

    py::class_<LinearExpansion<SpinResolvedSelectedONVBasis>>(module,
                                                              ("LinearExpansion_" + suffix).c_str(),
                                                              description.c_str())

        // CONSTRUCTORS

        .def_static(
            "Constant",
            [](const SpinResolvedSelectedONVBasis& onv_basis) {
                return LinearExpansion<SpinResolvedSelectedONVBasis>::Constant(onv_basis);
            },
            py::arg("onv_basis"),
            "Return a linear expansion with a normalized coefficient vector (i.e. all the coefficients are equal).")

        .def_static(
            "HartreeFock",
            [](const SpinResolvedSelectedONVBasis& onv_basis) {
                return LinearExpansion<SpinResolvedSelectedONVBasis>::HartreeFock(onv_basis);
            },
            py::arg("onv_basis"),
            "Return a linear expansion that represents the Hartree-Fock wave function.")

        .def_static(
            "Random",
            [](const SpinResolvedSelectedONVBasis& onv_basis) {
                return LinearExpansion<SpinResolvedSelectedONVBasis>::Random(onv_basis);
            },
            py::arg("onv_basis"),
            "Return a linear expansion with a random, normalized coefficient vector, with coefficients uniformly distributed in [-1, +1] before any normalization.")


        // PUBLIC METHODS

        .def(
            "calculate1DM",
            [](const LinearExpansion<SpinResolvedSelectedONVBasis>& linear_expansion) {
                return linear_expansion.calculate1DM();
            },
            "Return the one-electron density matrix (1-DM) for a selected spin-resolved wave function expansion.")

        .def(
            "calculate2DM",
            [](const LinearExpansion<SpinResolvedSelectedONVBasis>& linear_expansion) {
                return asNumpyArray(linear_expansion.calculate2DM().Eigen());
            },
            "Return the two-electron density matrix (2-DM) for a selected spin-resolved wave function expansion.")

        .def("coefficients",
             &LinearExpansion<SpinResolvedSelectedONVBasis>::coefficients,
             "Return the expansion coefficients of this linear expansion wave function model.")

        .def(
            "calculateSpinResolved1DM",
            [](const LinearExpansion<SpinResolvedSelectedONVBasis>& linear_expansion) {
                return linear_expansion.calculateSpinResolved1DM();
            },
            "Return the spin-resolved 1-DM.")

        .def(
            "calculateSpinResolved2DM",
            [](const LinearExpansion<SpinResolvedSelectedONVBasis>& linear_expansion) {
                return linear_expansion.calculateSpinResolved2DM();
            },
            "Return the spin resolved 2-DM.");
}


/**
 *  A template specialization for the binding of LinearExpansion<SpinUnresolvedONVBasis>, because it has an additional function to be bound.
 * 
 *  @param module               the Pybind11 module
 *  @param suffix               the suffix for the gqcpy class name, i.e. "LinearExpansion" + suffix
 *  @param description          the description for the gqcpy class
 */
template <>
void bindLinearExpansion<SpinUnresolvedONVBasis>(py::module& module, const std::string& suffix, const std::string& description) {

    py::class_<LinearExpansion<SpinUnresolvedONVBasis>>(module,
                                                        ("LinearExpansion_" + suffix).c_str(),
                                                        description.c_str())

        // CONSTRUCTORS

        .def_static(
            "Constant",
            [](const SpinUnresolvedONVBasis& onv_basis) {
                return LinearExpansion<SpinUnresolvedONVBasis>::Constant(onv_basis);
            },
            py::arg("onv_basis"),
            "Return a linear expansion with a normalized coefficient vector (i.e. all the coefficients are equal).")

        .def_static(
            "FromONVProjection",
            [](const SpinUnresolvedONV& onv_of, const GSpinorBasis<double, GTOShell>& spinor_basis_on, const GSpinorBasis<double, GTOShell>& spinor_basis_of) {
                return LinearExpansion<SpinUnresolvedONVBasis>::FromONVProjection(onv_of, spinor_basis_on, spinor_basis_of);
            },
            py::arg("onv_of"),
            py::arg("spinor_basis_on"),
            py::arg("spinor_basis_of"),
            "Create the linear expansion of the given spin-unresolved ONV that is expressed in the given GSpinorBasis, by projection onto the spin-resolved ONVs expressed with respect to another given GSpinorBasis.")

        .def_static(
            "HartreeFock",
            [](const SpinUnresolvedONVBasis& onv_basis) {
                return LinearExpansion<SpinUnresolvedONVBasis>::HartreeFock(onv_basis);
            },
            py::arg("onv_basis"),
            "Return a linear expansion that represents the Hartree-Fock wave function.")

        .def_static(
            "Random",
            [](const SpinUnresolvedONVBasis& onv_basis) {
                return LinearExpansion<SpinUnresolvedONVBasis>::Random(onv_basis);
            },
            py::arg("onv_basis"),
            "Return a linear expansion with a random, normalized coefficient vector, with coefficients uniformly distributed in [-1, +1] before any normalization.")


        // PUBLIC METHODS

        .def("coefficients",
             &LinearExpansion<SpinUnresolvedONVBasis>::coefficients,
             "Return the expansion coefficients of this linear expansion wave function model.");
}


void bindLinearExpansions(py::module& module) {

    bindLinearExpansion<SeniorityZeroONVBasis>(module, "SeniorityZero", "The linear expansion (configuration interaction) wave function model in a seniority-zero ONV basis.");
    bindLinearExpansion<SpinResolvedONVBasis>(module, "SpinResolved", "The linear expansion (configuration interaction) wave function model in a spin-resolved ONV basis.");
    bindLinearExpansion<SpinResolvedSelectedONVBasis>(module, "SpinResolvedSelected", "The linear expansion (configuration interaction) wave function model in a spin-resolved selected ONV basis.");
    bindLinearExpansion<SpinUnresolvedONVBasis>(module, "SpinUnresolved", "The linear expansion (configuration interaction) wave function model in a spin-unresolved ONV basis.");
}


}  // namespace gqcpy
