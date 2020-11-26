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

#include "QCModel/CC/CCSD.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


void bindQCModelCCSD(py::module& module) {

    py::class_<QCModel::CCSD<double>>(module, "QCModel_CCSD", "The CCSD wave function model.")

        // PUBLIC METHODS
        .def(
            "t1Amplitudes",
            &QCModel::CCSD<double>::t1Amplitudes,
            "Return these CCSD model parameters' T1-amplitudes")

        .def(
            "t2Amplitudes",
            &QCModel::CCSD<double>::t2Amplitudes,
            "Return these CCSD model parameters' T2-amplitudes");
}


}  // namespace gqcpy
