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
#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Basis/SpinorBasis/RSpinorBasis.hpp"
#include "Molecule/Molecule.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace py = pybind11;


namespace gqcpy {


void bindSQHamiltonian(py::module& module) {
    py::class_<GQCP::SQHamiltonian<double>>(module, "SQHamiltonian", "A class that represents a real, second-quantized Hamiltonian.")

        .def("Molecular", [] (const GQCP::RSpinorBasis<double, GQCP::GTOShell>& spinor_basis, const GQCP::Molecule& molecule) {
                return GQCP::SQHamiltonian<double>::Molecular(spinor_basis, molecule);
            },
            "Construct the molecular Hamiltonian in a given spinor basis."
        )

        .def("core", &GQCP::SQHamiltonian<double>::core, "Return the 'core' Hamiltonian, i.e. the total of the one-electron contributions to the Hamiltonian")

        .def("twoElectron", &GQCP::SQHamiltonian<double>::twoElectron, "Return the total of the two-electron contributions to the Hamiltonian")
    ;
}


}  // namespace gqcpy
