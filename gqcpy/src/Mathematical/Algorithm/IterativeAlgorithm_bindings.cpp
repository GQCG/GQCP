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
#include "Mathematical/Algorithm/IterativeAlgorithm.hpp"
#include "Mathematical/Optimization/Eigenproblem/EigenproblemEnvironment.hpp"
#include "Mathematical/Optimization/NonLinearEquation/NonLinearEquationEnvironment.hpp"
#include "QCMethod/CC/CCSDEnvironment.hpp"
#include "QCMethod/HF/GHF/GHFSCFEnvironment.hpp"
#include "QCMethod/HF/RHF/RHFSCFEnvironment.hpp"
#include "QCMethod/HF/UHF/UHFSCFEnvironment.hpp"
#include "Utilities/literals.hpp"

#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Since IterativeAlgorithm is a class template, we must provide bindings for each of its associated types. In order to avoid duplicating code, we use a templated binding approach.
 */

/**
 *  Bind an iterative algorithm to the given module.
 * 
 *  @tparam Environment         the type of the environment
 * 
 *  @param module               the Pybind11 module
 *  @param suffix               the suffix that the Python class should receive, i.e. "IterativeAlgorithm" + suffix
 *  @param description          the Python class description
 */
template <typename Environment>
void bindIterativeAlgorithm(py::module& module, const std::string& suffix, const std::string& description) {

    py::class_<IterativeAlgorithm<Environment>>(module,
                                                ("IterativeAlgorithm_" + suffix).c_str(),
                                                description.c_str())

        // PUBLIC METHODS

        .def(
            "description",
            &IterativeAlgorithm<Environment>::description,
            "Return a textual description of this iterative algorithm.")

        .def(
            "insert",
            [](IterativeAlgorithm<Environment>& algorithm, const FunctionalStep<Environment>& step, const size_t index) {
                algorithm.insert(step, index);
            },
            py::arg("step"),
            py::arg("index"),
            "Insert an algorithm step at the given index.")

        .def(
            "maximumNumberOfIterations",
            &IterativeAlgorithm<Environment>::maximumNumberOfIterations,
            "Return the maximum number of iterations the algorithm may perform")

        .def(
            "numberOfIterations",
            &IterativeAlgorithm<Environment>::numberOfIterations,
            "Return the number of iterations that have been performed")

        .def(
            "perform",
            [](IterativeAlgorithm<Environment>& algorithm, Environment& environment) {
                algorithm.perform(environment);
            },
            py::arg("environment"))


        .def(
            "remove",
            [](IterativeAlgorithm<Environment>& algorithm, const size_t index) {
                algorithm.remove(index);
            },
            py::arg("index"),
            "Remove the algorithm step at the given index.")

        .def(
            "replace",
            [](IterativeAlgorithm<Environment>& algorithm, const FunctionalStep<Environment>& step, const size_t index) {
                algorithm.replace(step, index);
            },
            py::arg("step"),
            py::arg("index"),
            "Replace an algorithm step at the given index.");
}


void bindIterativeAlgorithms(py::module& module) {

    bindIterativeAlgorithm<EigenproblemEnvironment>(module, "EigenproblemEnvironment", "An algorithm that performs iterations using an EigenproblemEnvironment.");
    bindIterativeAlgorithm<NonLinearEquationEnvironment<double>>(module, "NonLinearEquationEnvironment", "An algorithm that performs iterations using a NonLinearEquationEnvironment.");

    bindIterativeAlgorithm<RHFSCFEnvironment<double>>(module, "RHFSCFEnvironment", "An algorithm that performs iterations using an RHFSCFEnvironment.");
    bindIterativeAlgorithm<RHFSCFEnvironment<complex>>(module, "RHFSCFEnvironment_cd", "An algorithm that performs iterations using a complex RHFSCFEnvironment.");

    bindIterativeAlgorithm<UHFSCFEnvironment<double>>(module, "UHFSCFEnvironment", "An algorithm that performs iterations using an UHFSCFEnvironment.");
    bindIterativeAlgorithm<UHFSCFEnvironment<complex>>(module, "UHFSCFEnvironment_cd", "An algorithm that performs iterations using a complex UHFSCFEnvironment.");

    bindIterativeAlgorithm<GHFSCFEnvironment<double>>(module, "GHFSCFEnvironment_d", "An algorithm that performs iterations using a real GHFSCFEnvironment.");
    bindIterativeAlgorithm<GHFSCFEnvironment<complex>>(module, "GHFSCFEnvironment_cd", "An algorithm that performs iterations using a complex GHFSCFEnvironment.");

    bindIterativeAlgorithm<CCSDEnvironment<double>>(module, "CCSDEnvironment", "An algorithm that performs iterations using a CCSDEnvironment.");
}


}  // namespace gqcpy
