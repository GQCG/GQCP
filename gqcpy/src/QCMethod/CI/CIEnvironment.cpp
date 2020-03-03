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
#include "QCMethod/CI/CIEnvironment.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace py = pybind11;


namespace gqcpy {


void bindCIEnvironment(py::module& module) {

    auto module_ci_environment = module.def_submodule("CIEnvironment");

    module_ci_environment.def("Dense",
        [ ] (const GQCP::SQHamiltonian<double>& sq_hamiltonian, const GQCP::SeniorityZeroONVBasis& onv_basis) {
            return GQCP::CIEnvironment::Dense(sq_hamiltonian, onv_basis);
        },
        py::arg("sq_hamiltonian"),
        py::arg("onv_basis"),
        "Return an environment suitable for solving DOCI eigenvalue problems."
    );

    module_ci_environment.def("Dense",
        [ ] (const GQCP::HubbardHamiltonian<double>& hubbard_hamiltonian, const GQCP::SpinResolvedONVBasis& onv_basis) {
            return GQCP::CIEnvironment::Dense(hubbard_hamiltonian, onv_basis);
        },
        py::arg("hubbard_hamiltonian"),
        py::arg("onv_basis"),
        "Return an environment suitable for solving Hubbard-related eigenvalue problems."
    );

    module_ci_environment.def("Dense",
        [ ] (const GQCP::SQHamiltonian<double>& sq_hamiltonian, const GQCP::SpinResolvedONVBasis& onv_basis) {
            return GQCP::CIEnvironment::Dense(sq_hamiltonian, onv_basis);
        },
        py::arg("sq_hamiltonian"),
        py::arg("onv_basis"),
        "Return an environment suitable for solving spin-resolved FCI eigenvalue problems."
    );


    module_ci_environment.def("Iterative",
        [ ] (const GQCP::SQHamiltonian<double>& sq_hamiltonian, const GQCP::SeniorityZeroONVBasis& onv_basis, const GQCP::MatrixX<double>& V) {
            return GQCP::CIEnvironment::Iterative(sq_hamiltonian, onv_basis, V);
        },
        py::arg("sq_hamiltonian"),
        py::arg("onv_basis"),
        py::arg("V"),
        "Return an environment suitable for solving DOCI eigenvalue problems."
    );

    module_ci_environment.def("Iterative",
        [ ] (const GQCP::HubbardHamiltonian<double>& hubbard_hamiltonian, const GQCP::SpinResolvedONVBasis& onv_basis, const GQCP::MatrixX<double>& V) {
            return GQCP::CIEnvironment::Iterative(hubbard_hamiltonian, onv_basis, V);
        },
        py::arg("hubbard_hamiltonian"),
        py::arg("onv_basis"),
        py::arg("V"),
        "Return an environment suitable for solving Hubbard-related eigenvalue problems."
    );

    module_ci_environment.def("Iterative",
        [ ] (const GQCP::SQHamiltonian<double>& sq_hamiltonian, const GQCP::SpinResolvedONVBasis& onv_basis, const GQCP::MatrixX<double>& V) {
            return GQCP::CIEnvironment::Iterative(sq_hamiltonian, onv_basis, V);
        },
        py::arg("sq_hamiltonian"),
        py::arg("onv_basis"),
        py::arg("V"),
        "Return an environment suitable for solving spin-resolved FCI eigenvalue problems."
    );
}


}  // namespace gqcpy
