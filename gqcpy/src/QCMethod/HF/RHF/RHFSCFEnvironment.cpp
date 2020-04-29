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

#include "QCMethod/HF/RHF/RHFSCFEnvironment.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace py = pybind11;


namespace gqcpy {


void bindRHFSCFEnvironment(py::module& module) {
    py::class_<GQCP::RHFSCFEnvironment<double>>(module, "RHFSCFEnvironment", "An algorithm environment that can be used with standard RHF SCF solvers.")

        .def_static(
            "WithCoreGuess",
            [](const size_t N, const GQCP::SQHamiltonian<double>& sq_hamiltonian, const Eigen::MatrixXd& S) {  // use an itermediary Eigen matrix for the Python binding, since Pybind11 doesn't accept our types that are derived from Eigen::Matrix
                return GQCP::RHFSCFEnvironment<double>::WithCoreGuess(N, sq_hamiltonian, GQCP::QCMatrix<double> {S});
            },
            "Initialize an RHF SCF environment with an initial coefficient matrix that is obtained by diagonalizing the core Hamiltonian matrix.")


        .def_readwrite("N", &GQCP::RHFSCFEnvironment<double>::N);
}


}  // namespace gqcpy
