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

#pragma once

#include <unsupported/Eigen/CXX11/Tensor>

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Convert a rank-four Eigen::Tensor to a NumPy array.
 * 
 *  @param tensor       The `Eigen::Tensor` that should be converted to the NumPy array.
 * 
 *  @return The corresponding NumPy array.
 */
template <typename T>
py::array_t<T> asNumpyArray(const Eigen::Tensor<T, 4>& tensor) {

    // Implementation adapted from https://github.com/pybind/pybind11/issues/1377.
    const auto shape = tensor.dimensions();

    return py::array_t<T>(shape,
                          {shape[0] * shape[1] * shape[2] * sizeof(T), shape[1] * shape[2] * sizeof(T), shape[2] * sizeof(T), sizeof(T)},  // strides
                          tensor.data()                                                                                                    // data pointer
    );
}


}  // namespace gqcpy
