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

        // Bind the named constructor.
        .def_static(
            "WithCoreGuess",
            [](const size_t N, const GQCP::SQHamiltonian<double>& sq_hamiltonian, const Eigen::MatrixXd& S) {  // use an itermediary Eigen matrix for the Python binding, since Pybind11 doesn't accept our types that are derived from Eigen::Matrix
                return GQCP::RHFSCFEnvironment<double>::WithCoreGuess(N, sq_hamiltonian, GQCP::QCMatrix<double> {S});
            },
            "Initialize an RHF SCF environment with an initial coefficient matrix that is obtained by diagonalizing the core Hamiltonian matrix.")


        // Bind read-write members/properties, exposing intermediary environment variables to the Python interface.
        .def_readwrite("N", &GQCP::RHFSCFEnvironment<double>::N)

        .def_readwrite("orbital_energies", &GQCP::RHFSCFEnvironment<double>::orbital_energies)


        // The following properties require an explicit getter and setter because the setter can only accept native Eigen types.
        .def_property(
            "S",
            [](const GQCP::RHFSCFEnvironment<double>& environment) {
                return environment.S;
            },
            [](GQCP::RHFSCFEnvironment<double>& environment, const Eigen::MatrixXd& S) {
                environment.S = GQCP::QCMatrix<double>(S);
            })

        .def_property(
            "coefficient_matrices",
            [](const GQCP::RHFSCFEnvironment<double>& environment) {
                return environment.coefficient_matrices;
            },
            [](GQCP::RHFSCFEnvironment<double>& environment, const std::vector<Eigen::MatrixXd>& new_coefficient_matrices) {
                const std::deque<GQCP::TransformationMatrix<double>> coefficient_matrices {new_coefficient_matrices.begin(), new_coefficient_matrices.end()};
                environment.coefficient_matrices = coefficient_matrices;
            })

        .def_property(
            "density_matrices",
            [](const GQCP::RHFSCFEnvironment<double>& environment) {
                return environment.density_matrices;
            },
            [](GQCP::RHFSCFEnvironment<double>& environment, const std::vector<Eigen::MatrixXd>& new_density_matrices) {
                const std::deque<GQCP::OneRDM<double>> density_matrices {new_density_matrices.begin(), new_density_matrices.end()};
                environment.density_matrices = density_matrices;
            })

        .def_property(
            "fock_matrices",
            [](const GQCP::RHFSCFEnvironment<double>& environment) {
                return environment.fock_matrices;
            },
            [](GQCP::RHFSCFEnvironment<double>& environment, const std::vector<Eigen::MatrixXd>& new_fock_matrices) {
                const std::deque<GQCP::QCMatrix<double>> fock_matrices {new_fock_matrices.begin(), new_fock_matrices.end()};
                environment.fock_matrices = fock_matrices;
            })

        .def_property(
            "error_vectors",
            [](const GQCP::RHFSCFEnvironment<double>& environment) {
                return environment.error_vectors;
            },
            [](GQCP::RHFSCFEnvironment<double>& environment, const std::vector<Eigen::MatrixXd>& new_error_vectors) {
                const std::deque<GQCP::VectorX<double>> error_vectors {new_error_vectors.begin(), new_error_vectors.end()};
                environment.error_vectors = error_vectors;
            });
}


}  // namespace gqcpy
