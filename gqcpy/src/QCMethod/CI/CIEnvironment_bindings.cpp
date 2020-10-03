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

#include "QCMethod/CI/CIEnvironment.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace py = pybind11;


namespace gqcpy {


/**
 *  Bind a CI environment to a gqcpy submodule module.
 * 
 *  @tparam Hamiltonian             the type of the Hamiltonian
 *  @tparam ONVBasis                the type of the ONV basis
 * 
 *  @param submodule                the gqcpy.CIEnvironment submodule
 *  @param documentation            the documentation that should appear for the function
 */
template <typename Hamiltonian, typename ONVBasis>
void bindCIEnvironment(py::module& submodule, const std::string& documentation) {

    submodule.def(
        "Dense",
        [](const Hamiltonian& hamiltonian, const ONVBasis& onv_basis) {
            return GQCP::CIEnvironment::Dense(hamiltonian, onv_basis);
        },
        py::arg("hamiltonian"),
        py::arg("onv_basis"),
        documentation.c_str());

    submodule.def(
        "Iterative",
        [](const Hamiltonian& hamiltonian, const ONVBasis& onv_basis, const GQCP::MatrixX<double>& V) {
            return GQCP::CIEnvironment::Iterative(hamiltonian, onv_basis, V);
        },
        py::arg("hamiltonian"),
        py::arg("onv_basis"),
        py::arg("V"),
        documentation.c_str());
}


void bindCIEnvironments(py::module& module) {

    auto submodule = module.def_submodule("CIEnvironment");

    bindCIEnvironment<GQCP::RSQHamiltonian<double>, GQCP::SeniorityZeroONVBasis>(submodule, "Return an environment suitable for solving DOCI eigenvalue problems.");
    bindCIEnvironment<GQCP::HubbardHamiltonian<double>, GQCP::SpinResolvedONVBasis>(submodule, "Return an environment suitable for solving Hubbard-related eigenvalue problems.");
    bindCIEnvironment<GQCP::RSQHamiltonian<double>, GQCP::SpinResolvedONVBasis>(submodule, "Return an environment suitable for solving spin-resolved FCI eigenvalue problems.");
}


}  // namespace gqcpy
