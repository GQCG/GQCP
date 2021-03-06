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

#include "QCMethod/CC/CCSDSolver.hpp"

#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


void bindCCSDSolver(py::module& module) {
    py::class_<CCSDSolver<double>>(module, "CCSDSolver", "A factory class that can construct CCSD solvers in an easy way.")

        .def_static(
            "Plain",
            [](const double threshold, const size_t maximum_number_of_iterations) {
                return CCSDSolver<double>::Plain(threshold, maximum_number_of_iterations);
            },
            py::arg("threshold") = 1.0e-08,
            py::arg("maximum_number_of_iterations") = 128,
            "Return a plain CCSD solver that uses the norm of the difference of consecutive amplitudes as a convergence criterion.");
}


}  // namespace gqcpy
