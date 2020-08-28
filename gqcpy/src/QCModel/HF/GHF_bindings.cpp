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

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace py = pybind11;


namespace gqcpy {


void bindQCModelGHF(py::module& module) {
    py::class_<GQCP::QCModel::GHF<double>>(module, "QCModel_GHF", "The generalized Hartree-Fock wave function model.")

        // PUBLIC METHODS

        .def(
            "calculateOrthonormalBasis1DM",
            [](const GQCP::QCModel::GHF<double>& ghf_parameters) {
                return ghf_parameters.calculateOrthonormalBasis1DM();
            },
            "Return the 1-DM expressed in an orthonormal spinor basis related to these optimal GHF parameters.")

        .def(
            "calculateScalarBasis1DM",
            [](const GQCP::QCModel::GHF<double>& ghf_parameters) {
                return ghf_parameters.calculateScalarBasis1DM();
            },
            "Return the GHF 1-DM in the scalar/AO basis related to these optimal GHF parameters")

        .def("coefficientMatrix",
             &GQCP::QCModel::GHF<double>::coefficientMatrix,
             "Return the coefficient matrix that expresses every spatial orbital (as a column) in its underlying scalar basis.")

        .def("orbitalEnergies",
             &GQCP::QCModel::GHF<double>::orbitalEnergies,
             "Return the orbital energies.");
}


}  // namespace gqcpy
