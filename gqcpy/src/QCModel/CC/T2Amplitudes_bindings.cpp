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


namespace py = pybind11;


namespace gqcpy {

/**
 *  Convert a rank-four Eigen::Tensor as a numpy array
 * 
 *  @param tensor       the tensor that should be converted to the numpy array
 * 
 *  @return the corresponding numpy array
 */
template <typename T>
py::array_t<T> asNumpyArray(const Eigen::Tensor<T, 4>& tensor) {

    // Implementation adapted from https://github.com/pybind/pybind11/issues/1377
    const auto shape = tensor.dimensions();

    return py::array_t<T>(shape,
                          {shape[0] * shape[1] * shape[2] * sizeof(T), shape[1] * shape[2] * sizeof(T), shape[2] * sizeof(T), sizeof(T)},  // strides
                          tensor.data()                                                                                                    // data pointer
    );
}

void bindT2Amplitudes(py::module& module) {

    py::class_<GQCP::T2Amplitudes<double>>(module, "T2Amplitudes", "The coupled-cluster T2 amplitudes t_{ij}^{ab}.")

        // CONSTRUCTORS
        .def(py::init<>(
                 [](const Eigen::Tensor<double, 4>& T, const GQCP::OrbitalSpace& orbital_space) {
                     const auto t = orbital_space.createRepresentableObjectFor<double>(GQCP::OccupationType::k_occupied, GQCP::OccupationType::k_occupied, GQCP::OccupationType::k_virtual, GQCP::OccupationType::k_virtual, T);

                     return GQCP::T2Amplitudes<double>(t, orbital_space);
                 }),
             "Construct T2-amplitudes given their dense representation.")

        // PUBLIC METHODS
        .def(
            "asTensor",
            [](const GQCP::T2Amplitudes<double>& t2_amplitudes) {
                return asNumpyArray(t2_amplitudes.asImplicitRankFourTensorSlice().asTensor().Eigen());
            },
            "Return the T2-amplitudes as a NumPy array.")
        
        .def(
            "orbitalSpace",
            [](const GQCP::T2Amplitudes<double>& t2_amplitudes) {
                return t2_amplitudes.orbitalSpace();
            },
            "Return the OrbitalSpace associated with the T2-amplitudes.");
}


}  // namespace gqcpy
