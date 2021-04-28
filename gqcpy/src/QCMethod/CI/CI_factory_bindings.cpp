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

#include "ONVBasis/SeniorityZeroONVBasis.hpp"
#include "ONVBasis/SpinResolvedONVBasis.hpp"
#include "ONVBasis/SpinResolvedSelectedONVBasis.hpp"
#include "ONVBasis/SpinUnresolvedSelectedONVBasis.hpp"
#include "QCMethod/CI/CI.hpp"
#include "Utilities/complex.hpp"

#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  We want to provide an easier Python API for the construction of CI methods, e.g.
 *      gqcpy.CI(onv_basis)
 * 
 *  The following bindings implement this kind of factory behaviour.
 */

/**
 *  Bind a factory-like method for a CI method.
 * 
 *  @tparam Scalar              The scalar type for a CI expansion coefficient: real or complex.
 *  @tparam ONVBasis            The type of the ONV basis associated to the CI method.
 */
template <typename Scalar, typename ONVBasis>
void bindCIFactoryMethod(py::module& module) {

    module.def(
        "CI",
        [](const ONVBasis& onv_basis, const size_t number_of_states = 1) {
            return QCMethod::CI<Scalar, ONVBasis>(onv_basis, number_of_states);
        },
        "Return an appropriate CI method.",
        py::arg("onv_basis"),
        py::arg("number_of_states") = 1);
}


/**
 *  Bind all types of CI methods to the gqcpy module.
 */
void bindCIFactory(py::module& module) {

    // Bind real-valued CI methods.
    bindCIFactoryMethod<double, SeniorityZeroONVBasis>(module);
    bindCIFactoryMethod<double, SpinResolvedONVBasis>(module);
    bindCIFactoryMethod<double, SpinResolvedSelectedONVBasis>(module);

    bindCIFactoryMethod<double, SpinUnresolvedSelectedONVBasis>(module);


    // Bind complex-valued CI methods.
    module.def(
        "CI_cd",
        [](const SpinUnresolvedSelectedONVBasis& onv_basis, const size_t number_of_states = 1) {
            return QCMethod::CI<complex, SpinUnresolvedSelectedONVBasis>(onv_basis, number_of_states);
        },
        "Return an appropriate CI method.",
        py::arg("onv_basis"),
        py::arg("number_of_states") = 1);
}


}  // namespace gqcpy
