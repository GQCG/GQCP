// // This file is part of GQCG-GQCP.
// //
// // Copyright (C) 2017-2020  the GQCG developers
// //
// // GQCG-GQCP is free software: you can redistribute it and/or modify
// // it under the terms of the GNU Lesser General Public License as published by
// // the Free Software Foundation, either version 3 of the License, or
// // (at your option) any later version.
// //
// // GQCG-GQCP is distributed in the hope that it will be useful,
// // but WITHOUT ANY WARRANTY; without even the implied warranty of
// // MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// // GNU Lesser General Public License for more details.
// //
// // You should have received a copy of the GNU Lesser General Public License
// // along with GQCG-GQCP.  If not, see <http://www.gnu.org/licenses/>.

// #include "Mathematical/Algorithm/Algorithm.hpp"
// #include "Mathematical/Algorithm/IterativeAlgorithm.hpp"
// #include "ONVBasis/SeniorityZeroONVBasis.hpp"
// #include "ONVBasis/SpinResolvedONVBasis.hpp"
// #include "ONVBasis/SpinResolvedSelectedONVBasis.hpp"
// #include "QCMethod/CI/CI.hpp"
// #include "QCMethod/CI/CIEnvironment.hpp"

// #include <pybind11/pybind11.h>


// namespace py = pybind11;


// namespace gqcpy {


// /**
//  *  Since QCMethod::CI has a template argument for the representation of the ONVBasis, we'll have to bind each of them separately. In order to avoid duplicate code, we use a templated binding approach.
//  */

// /**
//  *  Bind a templated CI method.
//  *
//  *  @tparam ONVBasis            the scalar type of the SQOperator
//  *
//  *  @param module               the Pybind11 module
//  *  @param suffix               the suffix for the gqcpy class name, i.e. "SQOperator" + suffix
//  *  @param description          the description for the gqcpy class
//  */
// template <typename ONVBasis>
// void bindQCMethodCI(py::module& module, const std::string& suffix, const std::string& description) {
//     py::class_<GQCP::QCMethod::CI<ONVBasis>>(module,
//                                              ("CI" + suffix).c_str(),
//                                              description.c_str())

//         // CONSTRUCTORS

//         .def(py::init<const ONVBasis, const size_t>(),
//              py::arg("onv_basis"),
//              py::arg("number_of_states") = 1)


//         // PUBLIC METHODS

//         .def(
//             "optimize",
//             [](const GQCP::QCMethod::CI<ONVBasis>& qc_method, GQCP::Algorithm<GQCP::EigenproblemEnvironment>& solver, GQCP::EigenproblemEnvironment& environment) {
//                 return qc_method.optimize(solver, environment);
//             },
//             "Optimize the CI wave function model: find the linear expansion coefficients.")

//         .def(
//             "optimize",
//             [](const GQCP::QCMethod::CI<ONVBasis>& qc_method, GQCP::IterativeAlgorithm<GQCP::EigenproblemEnvironment>& solver, GQCP::EigenproblemEnvironment& environment) {
//                 return qc_method.optimize(solver, environment);
//             },
//             "Optimize the CI wave function model: find the linear expansion coefficients.");
// }


// void bindQCMethodCIs(py::module& module) {

//     bindQCMethodCI<GQCP::SeniorityZeroONVBasis>(module, "SeniorityZero", "Configuration interaction in a seniority-zero ONV basis.");
//     // bindQCMethodCI<GQCP::SpinResolvedFrozenONVBasis>(module, "SpinResolvedFrozen", "Configuration interaction in a frozen core spin-resolved ONV basis.");
//     bindQCMethodCI<GQCP::SpinResolvedONVBasis>(module, "SpinResolved", "Configuration interaction in a spin-resolved ONV basis.");
//     bindQCMethodCI<GQCP::SpinResolvedSelectedONVBasis>(module, "SpinResolvedSelected", "Configuration interaction in a spin-resolved selected ONV basis.");
// }


// }  // namespace gqcpy
