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

#include "Mathematical/Algorithm/FunctionalStep.hpp"
#include "Mathematical/Optimization/Eigenproblem/EigenproblemEnvironment.hpp"
#include "Mathematical/Optimization/NonLinearEquation/NonLinearEquationEnvironment.hpp"
#include "QCMethod/CC/CCSDEnvironment.hpp"
#include "QCMethod/HF/RHF/RHFSCFEnvironment.hpp"
#include "QCMethod/HF/UHF/UHFSCFEnvironment.hpp"

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Since FunctionalStep is a class template, we must provide bindings for each of its associated types. In order to avoid duplicating code, we use a templated binding approach.
 */

/**
 *  Bind a step to the given module.
 * 
 *  @tparam Environment         the type of the environment
 * 
 *  @param module               the Pybind11 module
 *  @param suffix               the suffix that the Python class should receive, i.e. "FunctionalStep" + suffix
 *  @param description          the Python class description
 */
template <typename Environment>
void bindFunctionalStep(py::module& module, const std::string& suffix, const std::string& description) {

    py::class_<FunctionalStep<Environment>>(module,
                                            ("FunctionalStep_" + suffix).c_str(),
                                            description.c_str())

        // CONSTRUCTORS

        .def(py::init<std::function<void(Environment&)>, const std::string&>(),
             py::arg("function"),
             py::arg("description") = "A custom functional step.");
}


void bindFunctionalSteps(py::module& module) {

    bindFunctionalStep<EigenproblemEnvironment>(module, "EigenproblemEnvironment", "A functional step that uses an EigenproblemEnvironment.");
    bindFunctionalStep<NonLinearEquationEnvironment<double>>(module, "NonLinearEquationEnvironment", "A functional step that uses a NonLinearEquationEnvironment.");

    bindFunctionalStep<RHFSCFEnvironment<double>>(module, "RHFSCFEnvironment", "A functional step that uses an RHFSCFEnvironment.");
    bindFunctionalStep<UHFSCFEnvironment<double>>(module, "UHFSCFEnvironment", "A functional step that uses an UHFSCFEnvironment.");

    bindFunctionalStep<CCSDEnvironment<double>>(module, "CCSDEnvironment", "A functional step that uses an CCSDEnvironment.");
}


}  // namespace gqcpy
