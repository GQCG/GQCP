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

#include "QCModel/HF/GHF.hpp"
#include "Utilities/aliases.hpp"
#include "gqcpy/include/interfaces.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


void bindQCModelGHF(py::module& module) {

    // Define Python class related to `real QCModel::GHF` and expose their interfaces.
    py::class_<QCModel::GHF<double>> py_QCModelGHF_d {module, "QCModel_GHF_d", "The generalized Hartree-Fock wave function model."};

    // Expose the `HartreeFock` interface.
    bindQCModelHartreeFockInterface(py_QCModelGHF_d);

    // Define Python class related to `real QCModel::GHF` and expose their interfaces.
    py::class_<QCModel::GHF<complex>> py_QCModelGHF_cd {module, "QCModel_GHF_cd", "The generalized Hartree-Fock wave function model."};

    // Expose the `HartreeFock` interface.
    bindQCModelHartreeFockInterface(py_QCModelGHF_cd);
}


}  // namespace gqcpy
