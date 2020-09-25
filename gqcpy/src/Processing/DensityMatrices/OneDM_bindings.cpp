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

#include "Processing/DensityMatrices/OneDM.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace py = pybind11;


namespace gqcpy {


void bindOneDM(py::module& module) {
    py::class_<GQCP::OneDM<double>>(module, "OneDM", "A single particle density matrix")

        // PUBLIC METHODS

        .def(
            "transformed",
            [](const Eigen::MatrixXd& D, const Eigen::MatrixXd& T) {
                return GQCP::OneDM<double> {D}.transformed(GQCP::TransformationMatrix<double> {T});
            },
            py::arg("T"),
            "Return the transformed density matrix.");
}


}  // namespace gqcpy