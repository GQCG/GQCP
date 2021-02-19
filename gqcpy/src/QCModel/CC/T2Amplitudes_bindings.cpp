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
#include "gqcpy/include/interfaces.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


void bindT2Amplitudes(py::module& module) {

    py::class_<T2Amplitudes<double>>(module, "T2Amplitudes", "The coupled-cluster T2 amplitudes t_{ij}^{ab}.")

        /*
         *  MARK: Access
         */
        .def(
            "orbitalSpace",
            &T2Amplitudes<double>::orbitalSpace,
            "Return the orbital space for these T2-amplitudes, which encapsulates the occupied-virtual separation.")

        .def(
            "get",
            [](const T2Amplitudes<double>& T2, const size_t i, const size_t j, const size_t a, const size_t b) {
                return T2(i, j, a, b);
            },
            py::arg("i"),
            py::arg("j"),
            py::arg("a"),
            py::arg("b"),
            "Return a read-only reference to the T2-amplitude t_{ij}^{ab}.")

        .def(
            "set",
            [](T2Amplitudes<double>& T2, const size_t i, const size_t j, const size_t a, const size_t b, const double value) {
                T2(i, j, a, b) = value;
            },
            py::arg("i"),
            py::arg("j"),
            py::arg("a"),
            py::arg("b"),
            py::arg("value"),
            "Change the T2-amplitude t_{ij}^{ab} to the given value.")


        /*
         *  MARK: Python-only methods
         */
        .def(
            "asTensor",
            [](const T2Amplitudes<double>& t2_amplitudes) {
                return asNumpyArray(t2_amplitudes.asImplicitRankFourTensorSlice().asTensor());
            },
            "Return the T2-amplitudes as a tensor.");
}


}  // namespace gqcpy
