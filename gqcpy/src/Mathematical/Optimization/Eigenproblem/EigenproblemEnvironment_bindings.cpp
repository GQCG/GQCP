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


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


void bindEigenproblemEnvironment(py::module& module) {

    py::class_<EigenproblemEnvironment>(module, "EigenproblemEnvironment", "An environment used to solve eigenvalue problems for self-adjoint matrices.")


        // Define read/write properties.
        .def_readwrite(
            "dimension",
            &EigenproblemEnvironment::dimension)


        // Bind properties with a custom setter (to allow for non-native Eigen types).
        .def_property(
            "A",
            [](const EigenproblemEnvironment& environment) {
                return environment.A;
            },
            [](EigenproblemEnvironment& environment, const Eigen::MatrixXd& new_A) {
                environment.A = SquareMatrix<double>(new_A);
            })

        .def_property(
            "diagonal",
            [](const EigenproblemEnvironment& environment) {
                return environment.diagonal;
            },
            [](EigenproblemEnvironment& environment, const Eigen::VectorXd& new_diagonal) {
                environment.diagonal = VectorX<double>(new_diagonal);
            })

        .def_property(
            "eigenvalues",
            [](const EigenproblemEnvironment& environment) {
                return environment.eigenvalues;
            },
            [](EigenproblemEnvironment& environment, const Eigen::VectorXd& new_eigenvalues) {
                environment.eigenvalues = VectorX<double>(new_eigenvalues);
            })

        .def_property(
            "eigenvectors",
            [](const EigenproblemEnvironment& environment) {
                return environment.eigenvectors;
            },
            [](EigenproblemEnvironment& environment, const Eigen::MatrixXd& new_eigenvectors) {
                environment.eigenvectors = MatrixX<double>(new_eigenvectors);
            })

        .def_property(
            "S",
            [](const EigenproblemEnvironment& environment) {
                return environment.S;
            },
            [](EigenproblemEnvironment& environment, const Eigen::MatrixXd& new_S) {
                environment.S = SquareMatrix<double>(new_S);
            })

        .def_property(
            "Lambda",
            [](const EigenproblemEnvironment& environment) {
                return environment.Lambda;
            },
            [](EigenproblemEnvironment& environment, const Eigen::VectorXd& new_Lambda) {
                environment.Lambda = VectorX<double>(new_Lambda);
            })

        .def_property(
            "Z",
            [](const EigenproblemEnvironment& environment) {
                return environment.Z;
            },
            [](EigenproblemEnvironment& environment, const Eigen::MatrixXd& new_Z) {
                environment.Z = MatrixX<double>(new_Z);
            })

        .def_property(
            "V",
            [](const EigenproblemEnvironment& environment) {
                return environment.V;
            },
            [](EigenproblemEnvironment& environment, const Eigen::MatrixXd& new_V) {
                environment.V = MatrixX<double>(new_V);
            })

        .def_property(
            "VA",
            [](const EigenproblemEnvironment& environment) {
                return environment.VA;
            },
            [](EigenproblemEnvironment& environment, const Eigen::MatrixXd& new_VA) {
                environment.VA = MatrixX<double>(new_VA);
            })

        .def_property(
            "X",
            [](const EigenproblemEnvironment& environment) {
                return environment.X;
            },
            [](EigenproblemEnvironment& environment, const Eigen::MatrixXd& new_X) {
                environment.X = MatrixX<double>(new_X);
            })

        .def_property(
            "R",
            [](const EigenproblemEnvironment& environment) {
                return environment.R;
            },
            [](EigenproblemEnvironment& environment, const Eigen::MatrixXd& new_R) {
                environment.R = MatrixX<double>(new_R);
            })

        .def_property(
            "Delta",
            [](const EigenproblemEnvironment& environment) {
                return environment.Delta;
            },
            [](EigenproblemEnvironment& environment, const Eigen::MatrixXd& new_Delta) {
                environment.Delta = MatrixX<double>(new_Delta);
            });
}


}  // namespace gqcpy
