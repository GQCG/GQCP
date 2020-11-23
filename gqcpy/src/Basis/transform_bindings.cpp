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

#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Basis/Transformations/RTransformation.hpp"
#include "Basis/Transformations/transform.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace py = pybind11;


namespace gqcpy {


void bindBasisTransform(py::module& module) {

    module.def(
        "transform",
        [](const Eigen::MatrixXd& T, GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell>& spinor_basis, GQCP::RSQHamiltonian<double>& sq_hamiltonian) {
            GQCP::transform(GQCP::RTransformation<double> {T}, spinor_basis, sq_hamiltonian);
        },
        py::arg("T"),
        py::arg("spinor_basis"),
        py::arg("sq_hamiltonian"));
}


}  // namespace gqcpy
