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

#include "Mathematical/Algorithm/Algorithm.hpp"
#include "Mathematical/Optimization/Eigenproblem/EigenproblemEnvironment.hpp"
#include "Mathematical/Optimization/LinearEquation/LinearEquationEnvironment.hpp"

#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Since Algorithm is a class template, we must provide bindings for each of its associated types. In order to avoid duplicating code, we use a templated binding approach.
 */

/**
 *  Bind an algorithm to the given module.
 * 
 *  @tparam Environment         the type of the environment
 * 
 *  @param module               the Pybind11 module
 *  @param suffix               the suffix that the Python class should receive, i.e. "Algorithm" + suffix
 *  @param description          the Python class description
 */
template <typename Environment>
void bindAlgorithm(py::module& module, const std::string& suffix, const std::string& description) {

    py::class_<Algorithm<Environment>>(module,
                                       ("Algorithm_" + suffix).c_str(),
                                       description.c_str())

        // PUBLIC METHODS

        .def("description",
             &Algorithm<Environment>::description,
             "Return a textual description of this algorithm.");
}


void bindAlgorithms(py::module& module) {

    bindAlgorithm<EigenproblemEnvironment>(module, "EigenproblemEnvironment", "An algorithm that only performs one collection of steps using an EigenproblemEnvironment.");
    bindAlgorithm<LinearEquationEnvironment<double>>(module, "LinearEquationEnvironment", "An algorithm that only performs one collection of steps using a LinearEquationEnvironment.");
}


}  // namespace gqcpy
