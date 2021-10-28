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
#include "QCMethod/NOCI/NOCI.hpp"
#include "Utilities/complex.hpp"

#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  We want to provide an easier Python API for the construction of NOCI methods, e.g.
 *      gqcpy.NOCI(non_orthogonal_basis)
 *
 *  The following bindings implement this kind of factory behaviour.
 */

/**
 *  Bind a factory-like method for a real NOCI method.
 *
 *  @tparam NonOrthogonalBasis            The type of the non-orthogonal basis associated to the NOCI method.
 */
template <typename NonOrthogonalBasis>
void bindRealNOCIFactoryMethod(py::module& module) {

    module.def(
        "NOCI_d",
        [](const NonOrthogonalBasis& non_orthogonal_basis, const size_t number_of_states = 1) {
            return QCMethod::NOCI<double, NonOrthogonalBasis>(non_orthogonal_basis, number_of_states);
        },
        "Return an appropriate NOCI method.",
        py::arg("non_orthogonal_basis"),
        py::arg("number_of_states") = 1);
}


/**
 *  Bind a factory-like method for a complex NOCI method.
 *
 *  @tparam NonOrthogonalBasis            The type of the non-orthogonal basis associated to the NOCI method.
 */
template <typename NonOrthogonalBasis>
void bindComplexNOCIFactoryMethod(py::module& module) {

    module.def(
        "NOCI_cd",
        [](const NonOrthogonalBasis& non_orthogonal_basis, const size_t number_of_states = 1) {
            return QCMethod::NOCI<complex, NonOrthogonalBasis>(non_orthogonal_basis, number_of_states);
        },
        "Return an appropriate NOCI method.",
        py::arg("non_orthogonal_basis"),
        py::arg("number_of_states") = 1);
}


/**
 *  Bind all types of NOCI methods to the gqcpy module.
 */
void bindNOCIFactory(py::module& module) {

    // Bind real-valued CI methods.
    bindRealNOCIFactoryMethod<GNonOrthogonalStateBasis<double>>(module);
    bindRealNOCIFactoryMethod<RNonOrthogonalStateBasis<double>>(module);
    bindRealNOCIFactoryMethod<UNonOrthogonalStateBasis<double>>(module);

    // Bind complex-valued CI methods.
    bindComplexNOCIFactoryMethod<GNonOrthogonalStateBasis<complex>>(module);
    bindComplexNOCIFactoryMethod<RNonOrthogonalStateBasis<complex>>(module);
    bindComplexNOCIFactoryMethod<UNonOrthogonalStateBasis<complex>>(module);
}


}  // namespace gqcpy
