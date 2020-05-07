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

#include "QCMethod/HF/UHF/UHFSCFEnvironment.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace py = pybind11;


namespace gqcpy {


void bindUHFSCFEnvironment(py::module& module) {
    py::class_<GQCP::UHFSCFEnvironment<double>>(module, "UHFSCFEnvironment", "An algorithm environment that can be used with standard UHF SCF solvers.")

        // Bind the constructor based on converged RHF model parameters.
        .def(py::init([](const GQCP::QCModel::RHF<double>& rhf_parameters, const GQCP::SQHamiltonian<double>& sq_hamiltonian, const Eigen::MatrixXd& S) {  // use an itermediary Eigen matrix for the Python binding, since Pybind11 doesn't accept our types that are derived from Eigen::Matrix
                 return GQCP::UHFSCFEnvironment<double>(rhf_parameters, sq_hamiltonian, GQCP::QCMatrix<double>(S));
             }),
             py::arg("rhf_parameters"),
             py::arg("sq_hamiltonian"),
             py::arg("S"),
             "A constructor that initializes the environment from converged RHF model parameters.")

        // Bind the named constructor.
        .def_static(
            "WithCoreGuess",
            [](const size_t N_alpha, const size_t N_beta, const GQCP::SQHamiltonian<double>& sq_hamiltonian, const Eigen::MatrixXd& S) {  // use an itermediary Eigen matrix for the Python binding, since Pybind11 doesn't accept our types that are derived from Eigen::Matrix
                return GQCP::UHFSCFEnvironment<double>::WithCoreGuess(N_alpha, N_beta, sq_hamiltonian, GQCP::QCMatrix<double>(S));
            },
            "Initialize an UHF SCF environment with initial coefficient matrices (equal for alpha and beta) that is obtained by diagonalizing the core Hamiltonian matrix.")


        // Bind read-write members/properties, exposing intermediary environment variables to the Python interface.
        .def_readwrite("N_alpha", &GQCP::UHFSCFEnvironment<double>::N_alpha)
        .def_readwrite("N_beta", &GQCP::UHFSCFEnvironment<double>::N_beta)

        .def_readwrite("orbital_energies_alpha", &GQCP::UHFSCFEnvironment<double>::orbital_energies_alpha)
        .def_readwrite("orbital_energies_beta", &GQCP::UHFSCFEnvironment<double>::orbital_energies_beta)


        // The following properties require an explicit getter and setter because the setter can only accept native Eigen types.
        .def_property(
            "S",
            [](const GQCP::UHFSCFEnvironment<double>& environment) {
                return environment.S;
            },
            [](GQCP::UHFSCFEnvironment<double>& environment, const Eigen::MatrixXd& S) {
                environment.S = GQCP::QCMatrix<double>(S);
            })


        .def_property(
            "coefficient_matrices_alpha",
            [](const GQCP::UHFSCFEnvironment<double>& environment) {
                return environment.coefficient_matrices_alpha;
            },
            [](GQCP::UHFSCFEnvironment<double>& environment, const std::vector<Eigen::MatrixXd>& new_coefficient_matrices_alpha) {
                const std::deque<GQCP::TransformationMatrix<double>> coefficient_matrices_alpha {new_coefficient_matrices_alpha.begin(), new_coefficient_matrices_alpha.end()};
                environment.coefficient_matrices_alpha = coefficient_matrices_alpha;
            })

        .def_property(
            "coefficient_matrices_beta",
            [](const GQCP::UHFSCFEnvironment<double>& environment) {
                return environment.coefficient_matrices_beta;
            },
            [](GQCP::UHFSCFEnvironment<double>& environment, const std::vector<Eigen::MatrixXd>& new_coefficient_matrices_beta) {
                const std::deque<GQCP::TransformationMatrix<double>> coefficient_matrices_beta {new_coefficient_matrices_beta.begin(), new_coefficient_matrices_beta.end()};
                environment.coefficient_matrices_beta = coefficient_matrices_beta;
            })

        .def_property(
            "density_matrices_alpha",
            [](const GQCP::UHFSCFEnvironment<double>& environment) {
                return environment.density_matrices_alpha;
            },
            [](GQCP::UHFSCFEnvironment<double>& environment, const std::vector<Eigen::MatrixXd>& new_density_matrices_alpha) {
                const std::deque<GQCP::OneRDM<double>> density_matrices_alpha {new_density_matrices_alpha.begin(), new_density_matrices_alpha.end()};
                environment.density_matrices_alpha = density_matrices_alpha;
            })

        .def_property(
            "density_matrices_beta",
            [](const GQCP::UHFSCFEnvironment<double>& environment) {
                return environment.density_matrices_beta;
            },
            [](GQCP::UHFSCFEnvironment<double>& environment, const std::vector<Eigen::MatrixXd>& new_density_matrices_beta) {
                const std::deque<GQCP::OneRDM<double>> density_matrices_beta {new_density_matrices_beta.begin(), new_density_matrices_beta.end()};
                environment.density_matrices_beta = density_matrices_beta;
            })


        .def("set_fock_matrices_alpha",
             [](GQCP::UHFSCFEnvironment<double>& environment, const std::vector<Eigen::MatrixXd>& new_fock_matrices_alpha) {
                 py::scoped_ostream_redirect stream {
                     std::cout,                                // std::ostream&
                     py::module::import("sys").attr("stdout")  // Python output
                 };

                 std::cout << "new_fock_matrices_alpha" << std::endl;
                 for (const auto& each : new_fock_matrices_alpha) {
                     std::cout << each << std::endl;
                 }
                 std::cout << std::endl;
                 //
             })


        .def("replace_current_fock_matrix_alpha",
             [](GQCP::UHFSCFEnvironment<double>& environment, const Eigen::MatrixXd& new_fock_matrix_alpha) {
                 environment.fock_matrices_alpha.pop_back();
                 environment.fock_matrices_alpha.push_back(GQCP::QCMatrix<double>(new_fock_matrix_alpha));
             })


        .def("replace_current_fock_matrix_beta",
             [](GQCP::UHFSCFEnvironment<double>& environment, const Eigen::MatrixXd& new_fock_matrix_beta) {
                 environment.fock_matrices_beta.pop_back();
                 environment.fock_matrices_beta.push_back(GQCP::QCMatrix<double>(new_fock_matrix_beta));
             })


        .def_property(
            "fock_matrices_alpha",
            [](const GQCP::UHFSCFEnvironment<double>& environment) {
                return environment.fock_matrices_alpha;
            },
            [](GQCP::UHFSCFEnvironment<double>& environment, const std::vector<Eigen::MatrixXd>& new_fock_matrices_alpha) {
                py::scoped_ostream_redirect stream {
                    std::cout,                                // std::ostream&
                    py::module::import("sys").attr("stdout")  // Python output
                };

                std::cout << "new_fock_matrices_alpha" << std::endl;
                for (const auto& each : new_fock_matrices_alpha) {
                    std::cout << each << std::endl;
                }
                std::cout << std::endl;
                //

                const std::deque<GQCP::QCMatrix<double>> fock_matrices_alpha {new_fock_matrices_alpha.begin(), new_fock_matrices_alpha.end()};
                std::cout << "fock_matrices_alpha" << std::endl;
                for (const auto& each : fock_matrices_alpha) {
                    std::cout << each << std::endl;
                }
                std::cout << std::endl;

                // environment.fock_matrices_alpha = fock_matrices_alpha;
            })

        .def_property(
            "fock_matrices_beta",
            [](const GQCP::UHFSCFEnvironment<double>& environment) {
                return environment.fock_matrices_beta;
            },
            [](GQCP::UHFSCFEnvironment<double>& environment, const std::vector<Eigen::MatrixXd>& new_fock_matrices_beta) {
                const std::deque<GQCP::QCMatrix<double>> fock_matrices_beta {new_fock_matrices_beta.begin(), new_fock_matrices_beta.end()};
                environment.fock_matrices_beta = fock_matrices_beta;
            })


        .def_property(
            "error_vectors_alpha",
            [](const GQCP::UHFSCFEnvironment<double>& environment) {
                return environment.error_vectors_alpha;
            },
            [](GQCP::UHFSCFEnvironment<double>& environment, const std::vector<Eigen::MatrixXd>& new_error_vectors_alpha) {
                const std::deque<GQCP::VectorX<double>> error_vectors_alpha {new_error_vectors_alpha.begin(), new_error_vectors_alpha.end()};
                environment.error_vectors_alpha = error_vectors_alpha;
            })

        .def_property(
            "error_vectors_beta",
            [](const GQCP::UHFSCFEnvironment<double>& environment) {
                return environment.error_vectors_beta;
            },
            [](GQCP::UHFSCFEnvironment<double>& environment, const std::vector<Eigen::MatrixXd>& new_error_vectors_beta) {
                const std::deque<GQCP::VectorX<double>> error_vectors_beta {new_error_vectors_beta.begin(), new_error_vectors_beta.end()};
                environment.error_vectors_beta = error_vectors_beta;
            });
}


}  // namespace gqcpy
