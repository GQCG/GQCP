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

#include "Basis/NonOrthogonalBasis/GNonOrthogonalStateBasis.hpp"
#include "Basis/NonOrthogonalBasis/RNonOrthogonalStateBasis.hpp"
#include "Basis/NonOrthogonalBasis/UNonOrthogonalStateBasis.hpp"
#include "Molecule/Molecule.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/NOCI/NOCIEnvironment.hpp"
#include "Utilities/complex.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Bind a real NOCI environment to a gqcpy submodule module.
 *
 *  @tparam Hamiltonian                       The type of the Hamiltonian.
 *  @tparam NonOrthogonalBasis                The type of the non-orthogonal basis.
 *
 *  @param submodule                The gqcpy.NOCIEnvironment submodule.
 *  @param documentation            The documentation that should appear for the function.
 */
template <typename Hamiltonian, typename NonOrthogonalBasis>
void bindRealNOCIEnvironment(py::module& submodule, const std::string& documentation) {

    submodule.def(
        "Dense_d",
        [](const Hamiltonian& hamiltonian, const NonOrthogonalBasis& non_orthogonal_basis, const Molecule& molecule) {
            return NOCIEnvironment::Dense(hamiltonian, non_orthogonal_basis, molecule);
        },
        py::arg("hamiltonian"),
        py::arg("non_orthogonal_basis"),
        py::arg("molecule"),
        documentation.c_str());
}


/**
 *  Bind a complex NOCI environment to a gqcpy submodule module.
 *
 *  @tparam Hamiltonian                       The type of the Hamiltonian.
 *  @tparam NonOrthogonalBasis                The type of the non-orthogonal basis.
 *
 *  @param submodule                The gqcpy.NOCIEnvironment submodule.
 *  @param documentation            The documentation that should appear for the function.
 */
template <typename Hamiltonian, typename NonOrthogonalBasis>
void bindComplexNOCIEnvironment(py::module& submodule, const std::string& documentation) {

    submodule.def(
        "Dense_cd",
        [](const Hamiltonian& hamiltonian, const NonOrthogonalBasis& non_orthogonal_basis, const Molecule& molecule) {
            return NOCIEnvironment::Dense(hamiltonian, non_orthogonal_basis, molecule);
        },
        py::arg("hamiltonian"),
        py::arg("non_orthogonal_basis"),
        py::arg("molecule"),
        documentation.c_str());
}


void bindNOCIEnvironments(py::module& module) {

    auto submodule = module.def_submodule("NOCIEnvironment");

    bindRealNOCIEnvironment<GSQHamiltonian<double>, GNonOrthogonalStateBasis<double>>(submodule, "Return an environment suitable for real-valued solving non-orthogonal configuration interaction in a basis that consists of 'generalized' states.");
    bindRealNOCIEnvironment<RSQHamiltonian<double>, RNonOrthogonalStateBasis<double>>(submodule, "Return an environment suitable for real-valued solving non-orthogonal configuration interaction in a basis that consists of 'generalized' states.");
    bindRealNOCIEnvironment<USQHamiltonian<double>, UNonOrthogonalStateBasis<double>>(submodule, "Return an environment suitable for real-valued solving non-orthogonal configuration interaction in a basis that consists of 'generalized' states.");

    bindComplexNOCIEnvironment<GSQHamiltonian<complex>, GNonOrthogonalStateBasis<complex>>(submodule, "Return an environment suitable for complex-valued solving non-orthogonal configuration interaction in a basis that consists of 'generalized' states.");
    bindComplexNOCIEnvironment<RSQHamiltonian<complex>, RNonOrthogonalStateBasis<complex>>(submodule, "Return an environment suitable for complex-valued solving non-orthogonal configuration interaction in a basis that consists of 'generalized' states.");
    bindComplexNOCIEnvironment<USQHamiltonian<complex>, UNonOrthogonalStateBasis<complex>>(submodule, "Return an environment suitable for complex-valued solving non-orthogonal configuration interaction in a basis that consists of 'generalized' states.");
}


}  // namespace gqcpy
