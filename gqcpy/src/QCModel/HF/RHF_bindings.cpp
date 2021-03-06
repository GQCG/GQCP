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

#include "QCModel/HF/RHF.hpp"
#include "Utilities/aliases.hpp"
#include "gqcpy/include/interfaces.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;

template <typename Scalar>
void bindQCModelRHF(py::module& module, const std::string& name, const std::string& description) {

    // Define Python classes related to `QCModel::RHF` and expose their interfaces.
    py::class_<QCModel::RHF<Scalar>> py_QCModelRHF {module, name.c_str(), description.c_str()};

    // Expose the actual Python bindings unique to RHF.
    py_QCModelRHF

        .def(
            "spinOrbitalEnergiesBlocked",
            &QCModel::RHF<Scalar>::spinOrbitalEnergiesBlocked,
            "Return all the spin-orbital energies, with the alpha spin-orbital energies appearing before the beta spin-orbital energies.")

        .def(
            "spinOrbitalEnergiesInterleaved",
            &QCModel::RHF<Scalar>::spinOrbitalEnergiesInterleaved,
            "Return all the spin-orbital energies, with the alpha spin-orbital energies appearing before the beta spin-orbital energies.");

    // Expose the `HartreeFock` interface.
    bindQCModelHartreeFockInterface(py_QCModelRHF);
}


/**
 *  Bind all types of `QCModel::RHF`s.
 */
void bindQCModelsRHF(py::module& module) {

    bindQCModelRHF<double>(module, "QCModel_RHF_d", "The real restricted Hartree-Fock wave function model.");
    bindQCModelRHF<complex>(module, "QCModel_RHF_cd", "The complex restricted Hartree-Fock wave function model.");
}


}  // namespace gqcpy
