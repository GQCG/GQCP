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
#include "ONVBasis/SeniorityZeroONVBasis.hpp"
#include "ONVBasis/SpinResolvedONVBasis.hpp"
#include "ONVBasis/SpinResolvedSelectedONVBasis.hpp"
#include "ONVBasis/SpinUnresolvedSelectedONVBasis.hpp"
#include "QCMethod/QCStructure.hpp"
#include "QCModel/CC/CCD.hpp"
#include "QCModel/CC/CCSD.hpp"
#include "QCModel/CI/LinearExpansion.hpp"
#include "QCModel/Geminals/AP1roG.hpp"
#include "QCModel/Geminals/vAP1roG.hpp"
#include "QCModel/HF/GHF.hpp"
#include "QCModel/HF/RHF.hpp"
#include "QCModel/HF/UHF.hpp"
#include "QCModel/NOCI/NOCIExpansion.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Since QCStructure is a class template, we must provide bindings for each of its associated types. In order to avoid duplicating code, we use a templated binding approach.
 */

/**
 *  Bind a quantum chemical model to the given module.
 *
 *  @tparam QCModel             the type of the quantum chemical model
 *
 *  @param module               the Pybind11 module
 *  @param suffix               the suffix that the Python class should receive, i.e. "QCStructure" + suffix
 *  @param description          the Python class description
 */
template <typename QCModel, typename Scalar = double>
void bindQCStructure(py::module& module, const std::string& suffix, const std::string& description) {

    py::class_<QCStructure<QCModel, Scalar>>(module,
                                             ("QCStructure_" + suffix).c_str(),
                                             description.c_str())

        // PUBLIC METHODS

        .def(
            "energy",
            [](const QCStructure<QCModel, Scalar>& qc_structure, const size_t i) {
                return qc_structure.energy(i);
            },
            py::arg("i") = 0,
            "Return the electronic energy corresponding to the i-th excited state.")

        .def(
            "groundStateEnergy",
            [](const QCStructure<QCModel, Scalar>& qc_structure) {
                return qc_structure.groundStateEnergy();
            },
            "Return the ground state electronic energy for this quantum chemical structure.")

        .def(
            "groundStateParameters",
            [](const QCStructure<QCModel, Scalar>& qc_structure) {
                return qc_structure.groundStateParameters();
            },
            "Return the ground state model parameters for this quantum chemical structure.")

        .def(
            "parameters",
            [](const QCStructure<QCModel, Scalar>& qc_structure, const size_t i) {
                return qc_structure.parameters(i);
            },
            py::arg("i") = 0,
            "Return the parameters corresponding to the i-th excited state.");
}


void bindQCStructures(py::module& module) {

    bindQCStructure<LinearExpansion<double, SeniorityZeroONVBasis>>(module, "LinearExpansion_SeniorityZero", "A quantum chemical structure for real-valued linear expansions in a seniority-zero ONV basis.");
    bindQCStructure<LinearExpansion<double, SpinResolvedONVBasis>>(module, "LinearExpansion_SpinResolved", "A quantum chemical structure for real-valued linear expansions in a spin-resolved ONV basis.");
    bindQCStructure<LinearExpansion<double, SpinResolvedSelectedONVBasis>>(module, "LinearExpansion_SpinResolvedSelected", "A quantum chemical structure for real-valued linear expansions in a spin-resolved selected ONV basis.");

    bindQCStructure<LinearExpansion<double, SpinUnresolvedSelectedONVBasis>>(module, "LinearExpansion_SpinUnresolvedSelected_d", "A quantum chemical structure for real-valued linear expansions in a spin-unresolved selected ONV basis.");
    bindQCStructure<LinearExpansion<complex, SpinUnresolvedSelectedONVBasis>, complex>(module, "LinearExpansion_SpinUnresolvedSelected_cd", "A quantum chemical structure for complex-valued linear expansions in a spin-unresolved selected ONV basis.");

    bindQCStructure<NOCIExpansion<double, GNonOrthogonalStateBasis<double>>>(module, "NOCI_GNonOrthogonalState_d", "A quantum chemical structure for real-valued expansions in a generalized non-orthogonal state basis.");
    bindQCStructure<NOCIExpansion<complex, GNonOrthogonalStateBasis<complex>>>(module, "NOCI_GNonOrthogonalState_cd", "A quantum chemical structure for complex-valued expansions in a generalized non-orthogonal state basis.");
    bindQCStructure<NOCIExpansion<double, RNonOrthogonalStateBasis<double>>>(module, "NOCI_RNonOrthogonalState_d", "A quantum chemical structure for real-valued expansions in a restricted non-orthogonal state basis.");
    bindQCStructure<NOCIExpansion<complex, RNonOrthogonalStateBasis<complex>>>(module, "NOCI_RNonOrthogonalState_cd", "A quantum chemical structure for complex-valued expansions in a restricted non-orthogonal state basis.");
    bindQCStructure<NOCIExpansion<double, UNonOrthogonalStateBasis<double>>>(module, "NOCI_UNonOrthogonalState_d", "A quantum chemical structure for real-valued expansions in a unrestricted non-orthogonal state basis.");
    bindQCStructure<NOCIExpansion<complex, UNonOrthogonalStateBasis<complex>>>(module, "NOCI_UNonOrthogonalState_cd", "A quantum chemical structure for complex-valued expansions in a unrestricted non-orthogonal state basis.");

    bindQCStructure<QCModel::AP1roG>(module, "AP1roG", "A quantum chemical structure for real-valued AP1roG parameters.");
    bindQCStructure<QCModel::vAP1roG>(module, "vAP1roG", "A quantum chemical structure for real-valued vAP1roG parameters.");

    bindQCStructure<QCModel::GHF<double>>(module, "GHF_d", "A quantum chemical structure for real-valued GHF parameters.");
    bindQCStructure<QCModel::GHF<complex>, complex>(module, "GHF_cd", "A quantum chemical structure for complex-valued GHF parameters.");

    bindQCStructure<QCModel::RHF<double>>(module, "RHF_d", "A quantum chemical structure for real-valued RHF parameters.");
    bindQCStructure<QCModel::RHF<complex>, complex>(module, "RHF_cd", "A quantum chemical structure for complex-valued RHF parameters.");

    bindQCStructure<QCModel::UHF<double>>(module, "UHF_d", "A quantum chemical structure for real-valued UHF parameters.");
    bindQCStructure<QCModel::UHF<complex>, complex>(module, "UHF_cd", "A quantum chemical structure for complex-valued UHF parameters.");

    bindQCStructure<QCModel::CCSD<double>>(module, "CCSD_d", "A quantum chemical structure for real-valued CCSD parameters.");
    bindQCStructure<QCModel::CCD<double>>(module, "CCD_d", "A quantum chemical structure for real-valued CCD parameters.");

    bindQCStructure<QCModel::CCSD<complex>, complex>(module, "CCSD_cd", "A quantum chemical structure for complex-valued CCSD parameters.");
    bindQCStructure<QCModel::CCD<complex>, complex>(module, "CCD_cd", "A quantum chemical structure for complex-valued CCD parameters.");
}


}  // namespace gqcpy
