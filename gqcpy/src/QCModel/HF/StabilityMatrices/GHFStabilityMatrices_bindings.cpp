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

#include "QCModel/HF/StabilityMatrices/GHFStabilityMatrices.hpp"
#include "Utilities/aliases.hpp"
#include "gqcpy/include/interfaces.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


void bindGHFStabilityMatrices(py::module& module) {

    // Define Python classes related to `QCModel::GHF` and expose their interfaces.
    py::class_<GHFStabilityMatrices<double>> py_GHFStabilityMatrices_d {module, "GHFStabilityMatrices_d", "The real generalized Hartree-Fock stability matrices."};

    bindQCModelHartreeFockStabilityInterface(py_GHFStabilityMatrices_d);


    // py::class_<GHFStabilityMatrices<complex>> py_GHFStabilityMatrices_cd {module, "GHFStabilityMatrices_cd", "The complex generalized Hartree-Fock stability matrices."};

    // bindQCModelHartreeFockStabilityInterface(py_GHFStabilityMatrices_cd);
}


}  // namespace gqcpy
