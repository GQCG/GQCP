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

#include "QCModel/HF/UHF.hpp"
#include "Utilities/aliases.hpp"
#include "gqcpy/include/interfaces.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;

template <typename Scalar>
void bindQCModelUHF(py::module& module, const std::string& name, const std::string& description) {

    // Define Python classes related to `QCModel::UHF` and expose their interfaces.
    py::class_<QCModel::UHF<Scalar>> py_QCModelUHF {module, name.c_str(), description.c_str()};

    // Expose the actual Python bindings unique to UHF.
    py_QCModelUHF

        .def(
            "numberOfSpinOrbitals",
            [](const QCModel::UHF<Scalar>& uhf_parameters, const Spin sigma) {
                return uhf_parameters.numberOfSpinOrbitals(sigma);
            },
            py::arg("sigma"),
            "Return the number of sigma spin-orbitals that these UHF model parameters describe.")

        .def(
            "numberOfSpinOrbitals",
            [](const QCModel::UHF<Scalar>& uhf_parameters) {
                return uhf_parameters.numberOfSpinOrbitals();
            },
            "Return the number of spin-orbitals that these UHF model parameters describe.")


        .def(
            "spinOrbitalEnergiesBlocked",
            &QCModel::UHF<Scalar>::spinOrbitalEnergiesBlocked,
            "Return all the spin-orbital energies, with the alpha spin-orbital energies appearing before the beta spin-orbital energies");

    // Expose the `HartreeFock` interface.
    bindQCModelHartreeFockInterface(py_QCModelUHF);
}

/**
 *  Bind all types of `QCModel::UHF`s.
 */
void bindQCModelsUHF(py::module& module) {

    bindQCModelUHF<double>(module, "QCModel_UHF_d", "The real unrestricted Hartree-Fock wave function model.");
    bindQCModelUHF<complex>(module, "QCModel_UHF_cd", "The complex unrestricted Hartree-Fock wave function model.");
}


}  // namespace gqcpy
