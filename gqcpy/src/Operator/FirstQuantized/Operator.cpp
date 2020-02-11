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
#include "Operator/FirstQuantized/Operator.hpp"

#include <pybind11/pybind11.h>


namespace py = pybind11;


namespace gqcpy {


void bindOperator(py::module& module) {

    py::class_<GQCP::NuclearRepulsionOperator>(module, "NuclearRepulsionOperator", "The nuclear repulsion operator.")

        .def("value",
            &GQCP::NuclearRepulsionOperator::value,
            "Return the scalar value of this nuclear repulsion operator."
        )
    ;


    py::class_<GQCP::Operator>(module, "Operator", "A class that is used to construct operators using static methods, much like a factory class.")

        .def_static("NuclearRepulsion",
            &GQCP::Operator::NuclearRepulsion,
            "Return a NuclearRepulsionOperator"
        )
    ;
}


}  // namespace gqcpy
