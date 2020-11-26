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

#include "QCModel/CC/T1Amplitudes.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


void bindT1Amplitudes(py::module& module) {

    py::class_<T1Amplitudes<double>>(module, "T1Amplitudes", "The coupled-cluster T1 amplitudes t_i^a.")

        // CONSTRUCTORS
        .def(py::init<>(
                 [](const Eigen::MatrixXd& T, const OrbitalSpace& orbital_space) {
                     const auto t = orbital_space.createRepresentableObjectFor<double>(OccupationType::k_occupied, OccupationType::k_virtual, T);

                     return T1Amplitudes<double>(t, orbital_space);
                 }),
             "Construct T1-amplitudes given their dense representation.")

        // PUBLIC METHODS
        .def(
            "asMatrix",
            [](const T1Amplitudes<double>& t1_amplitudes) {
                return t1_amplitudes.asImplicitMatrixSlice().asMatrix();
            },
            "Return the T1-amplitudes as a matrix.");
}


}  // namespace gqcpy
