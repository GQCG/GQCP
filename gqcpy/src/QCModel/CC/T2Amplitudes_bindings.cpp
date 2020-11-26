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

#include "QCModel/CC/T2Amplitudes.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


void bindT2Amplitudes(py::module& module) {

    py::class_<T2Amplitudes<double>>(module, "T2Amplitudes", "The coupled-cluster T2 amplitudes t_{ij}^{ab}.")

        // CONSTRUCTORS
        .def(py::init<>(
                 [](const Eigen::Tensor<double, 4>& T, const OrbitalSpace& orbital_space) {
                     const auto t = orbital_space.createRepresentableObjectFor<double>(OccupationType::k_occupied, OccupationType::k_occupied, OccupationType::k_virtual, OccupationType::k_virtual, T);

                     return T2Amplitudes<double>(t, orbital_space);
                 }),
             "Construct T2-amplitudes given their dense representation.")

        // PUBLIC METHODS
        .def(
            "asTensor",
            [](const T2Amplitudes<double>& t2_amplitudes) {
                return t2_amplitudes.asImplicitRankFourTensorSlice().asTensor();
            },
            "Return the T2-amplitudes as a tensor.");
}


}  // namespace gqcpy
