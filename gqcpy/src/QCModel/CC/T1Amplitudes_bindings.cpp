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

        /*
         *  MARK: Access
         */
        .def(
            "orbitalSpace",
            &T1Amplitudes<double>::orbitalSpace,
            "Return the orbital space for these T1-amplitudes, which encapsulates the occupied-virtual separation.")

        .def(
            "get",
            [](const T1Amplitudes<double>& T1, const size_t i, const size_t a) {
                return T1(i, a);
            },
            py::arg("i"),
            py::arg("a"),
            "Return a read-only reference to the T1-amplitude t_i^a.")

        .def(
            "set",
            [](T1Amplitudes<double>& T1, const size_t i, const size_t a, const double value) {
                T1(i, a) = value;
            },
            py::arg("i"),
            py::arg("a"),
            py::arg("value"),
            "Change the T1-amplitude t_i^a to the given value.")


        /*
         *  MARK: Python-only methods
         */
        .def(
            "asMatrix",
            [](const T1Amplitudes<double>& t1_amplitudes) {
                return t1_amplitudes.asImplicitMatrixSlice().asMatrix();
            },
            "Return the T1-amplitudes as a matrix.");
}


}  // namespace gqcpy
