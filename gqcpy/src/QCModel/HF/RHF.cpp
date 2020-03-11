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
#include "QCModel/HF/RHF.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace py = pybind11;


namespace gqcpy {


void bindQCModelRHF(py::module& module) {
    py::class_<GQCP::QCModel::RHF<double>>(module, "QCModel_RHF", "The restricted Hartree-Fock wave function model.")

        .def("calculateOrthonormalBasis1RDM",
            [ ] (const GQCP::QCModel::RHF<double>& rhf_parameters) {
                return rhf_parameters.calculateOrthonormalBasis1RDM();
            },
            "Return the 1-RDM expressed in an orthonormal spinor basis related to these optimal RHF parameters."
        )

        .def("calculateScalarBasis1RDM",
            [ ] (const GQCP::QCModel::RHF<double>& rhf_parameters) {
                return rhf_parameters.calculateScalarBasis1RDM();
            },
            "Return the RHF 1-RDM in the scalar/AO basis related to these optimal RHF parameters"
        )

        .def("coefficientMatrix",
            &GQCP::QCModel::RHF<double>::coefficientMatrix,
            "Return the coefficient matrix that expresses every spatial orbital (as a column) in its underlying scalar basis."
        )

        .def("orbitalEnergies",
            &GQCP::QCModel::RHF<double>::orbitalEnergies,
            "Return the orbital energies."
        )
    ;
}



}  // namespace gqcpy
