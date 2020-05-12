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

#include "Mathematical/Optimization/Eigenproblem/EigenproblemEnvironment.hpp"

#include <pybind11/pybind11.h>


namespace py = pybind11;


namespace gqcpy {


void bindEigenproblemEnvironment(py::module& module) {

    py::class_<GQCP::EigenproblemEnvironment>(module, "EigenproblemEnvironment", "An environment used to solve eigenvalue problems for self-adjoint matrices.")


        // Define read/write properties.
        .def_readwrite(
            "dimension",
            &GQCP::EigenproblemEnvironment::dimension)


        // Bind properties with a custom setter (to allow for non-native Eigen types).
        .def_property(
            "A",
            [](const GQCP::EigenproblemEnvironment& environment) {
                return environment.A;
            },
            [](GQCP::EigenproblemEnvironment& environment, const Eigen::MatrixXd& new_A) {
                environment.A = GQCP::SquareMatrix<double>(new_A);
            })

        .def_property(
            "diagonal",
            [](const GQCP::EigenproblemEnvironment& environment) {
                return environment.diagonal;
            },
            [](GQCP::EigenproblemEnvironment& environment, const Eigen::VectorXd& new_diagonal) {
                environment.diagonal = GQCP::VectorX<double>(new_diagonal);
            })

        .def_property(
            "eigenvalues",
            [](const GQCP::EigenproblemEnvironment& environment) {
                return environment.eigenvalues;
            },
            [](GQCP::EigenproblemEnvironment& environment, const Eigen::VectorXd& new_eigenvalues) {
                environment.eigenvalues = GQCP::VectorX<double>(new_eigenvalues);
            })

        .def_property(
            "eigenvectors",
            [](const GQCP::EigenproblemEnvironment& environment) {
                return environment.eigenvectors;
            },
            [](GQCP::EigenproblemEnvironment& environment, const Eigen::MatrixXd& new_eigenvectors) {
                environment.eigenvectors = GQCP::MatrixX<double>(new_eigenvectors);
            })

        .def_property(
            "S",
            [](const GQCP::EigenproblemEnvironment& environment) {
                return environment.S;
            },
            [](GQCP::EigenproblemEnvironment& environment, const Eigen::MatrixXd& new_S) {
                environment.S = GQCP::SquareMatrix<double>(new_S);
            })

        .def_property(
            "Lambda",
            [](const GQCP::EigenproblemEnvironment& environment) {
                return environment.Lambda;
            },
            [](GQCP::EigenproblemEnvironment& environment, const Eigen::VectorXd& new_Lambda) {
                environment.Lambda = GQCP::VectorX<double>(new_Lambda);
            })

        .def_property(
            "Z",
            [](const GQCP::EigenproblemEnvironment& environment) {
                return environment.Z;
            },
            [](GQCP::EigenproblemEnvironment& environment, const Eigen::MatrixXd& new_Z) {
                environment.Z = GQCP::MatrixX<double>(new_Z);
            })

        .def_property(
            "V",
            [](const GQCP::EigenproblemEnvironment& environment) {
                return environment.V;
            },
            [](GQCP::EigenproblemEnvironment& environment, const Eigen::MatrixXd& new_V) {
                environment.V = GQCP::MatrixX<double>(new_V);
            })

        .def_property(
            "VA",
            [](const GQCP::EigenproblemEnvironment& environment) {
                return environment.VA;
            },
            [](GQCP::EigenproblemEnvironment& environment, const Eigen::MatrixXd& new_VA) {
                environment.VA = GQCP::MatrixX<double>(new_VA);
            })

        .def_property(
            "X",
            [](const GQCP::EigenproblemEnvironment& environment) {
                return environment.X;
            },
            [](GQCP::EigenproblemEnvironment& environment, const Eigen::MatrixXd& new_X) {
                environment.X = GQCP::MatrixX<double>(new_X);
            })

        .def_property(
            "R",
            [](const GQCP::EigenproblemEnvironment& environment) {
                return environment.R;
            },
            [](GQCP::EigenproblemEnvironment& environment, const Eigen::MatrixXd& new_R) {
                environment.R = GQCP::MatrixX<double>(new_R);
            })

        .def_property(
            "Delta",
            [](const GQCP::EigenproblemEnvironment& environment) {
                return environment.Delta;
            },
            [](GQCP::EigenproblemEnvironment& environment, const Eigen::MatrixXd& new_Delta) {
                environment.Delta = GQCP::MatrixX<double>(new_Delta);
            });
}


}  // namespace gqcpy
