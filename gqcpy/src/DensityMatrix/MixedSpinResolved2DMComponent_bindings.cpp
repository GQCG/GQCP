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

#include "DensityMatrix/MixedSpinResolved2DMComponent.hpp"
#include "gqcpy/include/utilities.hpp"

#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Register `MixedSpinResolved2DMComponent_d` to the gqcpy module and expose a part of its C++ interface to Python.
 * 
 *  @param module           The Pybind11 module in which `MixedSpinResolved2DMComponent_d` should be registered.
 */
void bindMixedSpinResolved2DMComponent(py::module& module) {

    // Define the Python class for `MixedSpinResolved2DMComponent`.
    py::class_<MixedSpinResolved2DMComponent<double>> py_MixedSpinResolved2DMComponent_d {module, "MixedSpinResolved2DMComponent_d", "One of the mixed (i.e. alpha-beta or beta-alpha) spin components of a spin-resolved 2-DM."};


    py_MixedSpinResolved2DMComponent_d
        .def(
            "asArray",
            [](const MixedSpinResolved2DMComponent<double>& d) {
                return asNumpyArray(d.Eigen());
            });
}


}  // namespace gqcpy
