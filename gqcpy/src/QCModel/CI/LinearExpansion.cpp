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
#include "ONVBasis/SeniorityZeroONVBasis.hpp"
#include "ONVBasis/SpinResolvedFrozenONVBasis.hpp"
#include "ONVBasis/SpinResolvedONVBasis.hpp"
#include "ONVBasis/SpinResolvedSelectedONVBasis.hpp"
#include "QCModel/CI/LinearExpansion.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace py = pybind11;


namespace gqcpy {


/**
 *  Since LinearExpansion has a template argument for the representation of the ONVBasis, we'll have to bind each of them separately. In order to avoid duplicate code, we use a templated binding approach.
 */

/**
 *  Bind a templated CI method.
 * 
 *  @tparam ONVBasis            the scalar type of the SQOperator
 * 
 *  @param module               the Pybind11 module
 *  @param suffix               the suffix for the gqcpy class name, i.e. "SQOperator" + suffix
 *  @param description          the description for the gqcpy class
 */
template <typename ONVBasis>
void bindLinearExpansion(py::module& module, const std::string& suffix, const std::string& description) {
    py::class_<GQCP::LinearExpansion<ONVBasis>>(module,
        ("LinearExpansion" + suffix).c_str(),
        description.c_str()
    )

        .def("coefficients",
            &GQCP::LinearExpansion<ONVBasis>::coefficients,
            "Return the expansion coefficients of this linear expansion wave function model."
        )
    ;
}


/**
 *  A template specialization for the binding of GQCP::LinearExpansion<GQCP::SeniorityZeroONVBasis>, because it has an additional function to be bound.
 * 
 *  @param module               the Pybind11 module
 *  @param suffix               the suffix for the gqcpy class name, i.e. "SQOperator" + suffix
 *  @param description          the description for the gqcpy class
 */
template <>
void bindLinearExpansion<GQCP::SeniorityZeroONVBasis>(py::module& module, const std::string& suffix, const std::string& description) {

    py::class_<GQCP::LinearExpansion<GQCP::SeniorityZeroONVBasis>>(module,
        ("LinearExpansion" + suffix).c_str(),
        description.c_str()
    )

        .def("calculate1DM",
            [ ] (const GQCP::LinearExpansion<GQCP::SeniorityZeroONVBasis>& linear_expansion) {
                return linear_expansion.calculate1DM();
            },
            "Return the one-electron density matrix (1-DM) for a seniority-zero wave function expansion."
        )

        .def("coefficients",
            &GQCP::LinearExpansion<GQCP::SeniorityZeroONVBasis>::coefficients,
            "Return the expansion coefficients of this linear expansion wave function model."
        )
    ;
}


void bindLinearExpansions(py::module& module) {

    bindLinearExpansion<GQCP::SeniorityZeroONVBasis>(module, "SeniorityZero", "The linear expansion (configuration interaction) wave function model in a seniority-zero ONV basis.");
    bindLinearExpansion<GQCP::SpinResolvedFrozenONVBasis>(module, "SpinResolvedFrozen", "The linear expansion (configuration interaction) wave function model in a frozen core spin-resolved ONV basis.");
    bindLinearExpansion<GQCP::SpinResolvedONVBasis>(module, "SpinResolved", "The linear expansion (configuration interaction) wave function model in a spin-resolved ONV basis.");
    bindLinearExpansion<GQCP::SpinResolvedSelectedONVBasis>(module, "SpinResolvedSelected", "The linear expansion (configuration interaction) wave function model in a spin-resolved selected ONV basis.");
}



}  // namespace gqcpy
