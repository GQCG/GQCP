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

#include "Mathematical/Algorithm/IterativeAlgorithm.hpp"
#include "QCMethod/HF/GHF/GHF.hpp"
#include "QCMethod/HF/GHF/GHFSCFSolver.hpp"
#include "Utilities/aliases.hpp"

#include <pybind11/pybind11.h>


namespace py = pybind11;


namespace gqcpy {


template <typename Scalar>
void bindQCMethodGHF(py::module& module, const std::string& suffix) {
    py::class_<GQCP::QCMethod::GHF<Scalar>>(module,
                                            ("GHF_" + suffix).c_str(),
                                            "The generalized Hartree-Fock quantum chemical method.")

        // PUBLIC METHODS

        .def_static(
            "optimize",
            [](GQCP::IterativeAlgorithm<GQCP::GHFSCFEnvironment<Scalar>>& solver, GQCP::GHFSCFEnvironment<Scalar>& environment) {
                return GQCP::QCMethod::GHF<Scalar>().optimize(solver, environment);
            },
            py::arg("solver"),
            py::arg("environment"),
            "Optimize the GHF wave function model.");
}


void bindQCMethodsGHF(py::module& module) {

    bindQCMethodGHF<double>(module, "d");          // suffix 'd' for the class name
    bindQCMethodGHF<GQCP::complex>(module, "cd");  // suffix 'cd' for the class name: 'complex double'
}


}  // namespace gqcpy
