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
#include "ONVBasis/SpinResolvedSelectedONVBasis.hpp"
#include "ONVBasis/SpinUnresolvedSelectedONVBasis.hpp"
#include "Operator/SecondQuantized/ModelHamiltonian/HubbardHamiltonian.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/CI/CIEnvironment.hpp"
#include "Utilities/complex.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Bind the full CI environment to a gqcpy submodule module.
 *
 *  @tparam Hamiltonian             the type of the Hamiltonian
 *  @tparam ONVBasis                the type of the ONV basis
 *
 *  @param submodule                the gqcpy.CIEnvironment submodule
 *  @param documentation            the documentation that should appear for the function
 */
template <typename Hamiltonian, typename ONVBasis>
void bindFullCIEnvironment(py::module& submodule, const std::string& documentation) {

    submodule.def(
        "Dense",
        [](const Hamiltonian& hamiltonian, const ONVBasis& onv_basis) {
            return CIEnvironment::Dense(hamiltonian, onv_basis);
        },
        py::arg("hamiltonian"),
        py::arg("onv_basis"),
        documentation.c_str());

    submodule.def(
        "Iterative",
        [](const Hamiltonian& hamiltonian, const ONVBasis& onv_basis, const MatrixX<double>& V) {
            return CIEnvironment::Iterative(hamiltonian, onv_basis, V);
        },
        py::arg("hamiltonian"),
        py::arg("onv_basis"),
        py::arg("V"),
        documentation.c_str());
}


/**
 *  Bind the Dense only CI environment to a gqcpy submodule module.
 *
 *  @tparam Hamiltonian             the type of the Hamiltonian
 *  @tparam ONVBasis                the type of the ONV basis
 *
 *  @param submodule                the gqcpy.CIEnvironment submodule
 *  @param documentation            the documentation that should appear for the function
 */
template <typename Hamiltonian, typename ONVBasis>
void bindDenseOnlyCIEnvironment(py::module& submodule, const std::string& documentation) {

    submodule.def(
        "Dense",
        [](const Hamiltonian& hamiltonian, const ONVBasis& onv_basis) {
            return CIEnvironment::Dense(hamiltonian, onv_basis);
        },
        py::arg("hamiltonian"),
        py::arg("onv_basis"),
        documentation.c_str());
}


void bindCIEnvironments(py::module& module) {

    auto submodule = module.def_submodule("CIEnvironment");

    bindFullCIEnvironment<RSQHamiltonian<double>, SpinResolvedONVBasis>(submodule, "Return an environment suitable for solving spin-resolved FCI eigenvalue problems.");
    bindFullCIEnvironment<USQHamiltonian<double>, SpinResolvedONVBasis>(submodule, "Return an environment suitable for solving spin-resolved FCI eigenvalue problems.");
    bindFullCIEnvironment<HubbardHamiltonian<double>, SpinResolvedONVBasis>(submodule, "Return an environment suitable for solving Hubbard problems.");

    bindFullCIEnvironment<RSQHamiltonian<double>, SpinResolvedSelectedONVBasis>(submodule, "Return an environment suitable for solving spin-resolved selected eigenvalue problems.");
    bindFullCIEnvironment<USQHamiltonian<double>, SpinResolvedSelectedONVBasis>(submodule, "Return an environment suitable for solving spin-resolved selected eigenvalue problems.");

    bindFullCIEnvironment<RSQHamiltonian<double>, SeniorityZeroONVBasis>(submodule, "Return an environment suitable for solving seniority zero eigenvalue problems.");

    // We haven't implemented matrix-vector products for `SpinUnresolved(Selected)ONVBasis`, so we can't use the `Iterative` APIs.
    bindDenseOnlyCIEnvironment<GSQHamiltonian<double>, SpinUnresolvedONVBasis>(submodule, "Return an environment suitable for solving spin-unresolved eigenvalue problems.");
    bindDenseOnlyCIEnvironment<GSQHamiltonian<double>, SpinUnresolvedSelectedONVBasis>(submodule, "Return an environment suitable for solving spin-unresolved selected eigenvalue problems.");

    submodule.def(
        "Dense_cd",
        [](const GSQHamiltonian<complex>& hamiltonian, const SpinUnresolvedSelectedONVBasis& onv_basis) {
            return CIEnvironment::Dense(hamiltonian, onv_basis);
        },
        py::arg("hamiltonian"),
        py::arg("onv_basis"),
        "Return an environment suitable for solving spin-unresolved selected eigenvalue problems.");
}


}  // namespace gqcpy
