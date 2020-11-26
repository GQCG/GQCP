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

#include "Processing/Properties/properties.hpp"
#include "Processing/Properties/vAP1roGElectricalResponseSolver.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


void bindvAP1roGElectricalResponseSolver(py::module& module) {

    py::class_<vAP1roGElectricalResponseSolver>(module, "vAP1roGElectricalResponseSolver", "A class whose instances can solve the response equations for vAP1roG.")

        // CONSTRUCTORS

        .def(py::init<const QCModel::vAP1roG>(),
             py::arg("vap1rog"))


        // PUBLIC METHODS

        .def(
            "calculateWaveFunctionResponse",
            [](const vAP1roGElectricalResponseSolver& response_solver, const RSQHamiltonian<double>& sq_hamiltonian, const VectorRSQOneElectronOperator<double>& dipole_op) {
                return response_solver.calculateWaveFunctionResponse(sq_hamiltonian, dipole_op);
            },
            "Solve the parameter-linear response equations and return the wave function response.",
            py::arg("sq_hamiltonian"),
            py::arg("dipole_op"))

        .def(
            "calculateMultiplierResponse",
            [](const vAP1roGElectricalResponseSolver& response_solver, const RSQHamiltonian<double>& sq_hamiltonian, const VectorRSQOneElectronOperator<double>& dipole_op, const Eigen::Matrix<double, Eigen::Dynamic, 3>& x) {
                return response_solver.calculateMultiplierResponse(sq_hamiltonian, dipole_op, x);
            },
            "Solve the multiplier-linear response equations and return the wave function response.")

        .def(
            "calculateElectricPolarizability",
            [](const vAP1roGElectricalResponseSolver& response_solver, const Eigen::Matrix<double, Eigen::Dynamic, 3>& x, const VectorRSQOneElectronOperator<double>& dipole_op, const Eigen::Matrix<double, Eigen::Dynamic, 3>& y) {
                const auto F_p = response_solver.calculateParameterResponseForce(dipole_op);
                const auto A_lambda = response_solver.calculateExplicitMultiplierResponseForce(dipole_op);

                return calculateElectricPolarizability(F_p, x, A_lambda, y);
            },
            "Calculate the AP1roG electric polarizability",
            py::arg("x"),
            py::arg("dipole_op"),
            py::arg("y"));
}


}  // namespace gqcpy
