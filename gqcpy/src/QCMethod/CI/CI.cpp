// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2019  the GQCG developers
// 
// GQCG-gqcp is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// GQCG-gqcp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-gqcp.  If not, see <http://www.gnu.org/licenses/>.
// 
#include "QCMethod/CI/CI.hpp"

#include "Mathematical/Algorithm/IterativeAlgorithm.hpp"
#include "ONVBasis/SpinResolvedONVBasis.hpp"
#include "QCMethod/CI/CIEnvironment.hpp"


#include <pybind11/pybind11.h>


namespace py = pybind11;


namespace gqcpy {


void bindQCMethodCI(py::module& module) {
    py::class_<GQCP::QCMethod::CI>(module, "CI", "The configuration interaction quantum chemical method.")

        .def(py::init<const GQCP::SpinResolvedONVBasis, const size_t>(),
            py::arg("onv_basis"),
            py::arg("number_of_states") = 1
        )

        .def("optimize",
            [] (const GQCP::QCMethod::CI& qc_method, GQCP::Algorithm<GQCP::EigenproblemEnvironment>& solver, GQCP::EigenproblemEnvironment& environment) {
                return qc_method.optimize(solver, environment);
            },
            "Optimize the CI wave function model: find the linear expansion coefficients."
        )

        .def("optimize",
            [] (const GQCP::QCMethod::CI& qc_method, GQCP::IterativeAlgorithm<GQCP::EigenproblemEnvironment>& solver, GQCP::EigenproblemEnvironment& environment) {
                return qc_method.optimize(solver, environment);
            },
            "Optimize the CI wave function model: find the linear expansion coefficients."
        )
    ;
}


}  // namespace gqcpy
