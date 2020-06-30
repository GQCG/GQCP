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

#include "QCMethod/CC/CCSDEnvironment.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace py = pybind11;


namespace gqcpy {


/**
 *  Convert a NumPy array to a rank-four Eigen::Tensor since the automatic conversion from a NumPy array to an Eigen::Tensor is not supported yet.
 * 
 *  @param array       the NumPy array that should be converted to a tensor
 * 
 *  @return the corresponding rank-four Eigen::Tensor
 */
template <typename T>
const Eigen::Tensor<T, 4> asTensor(const py::array_t<T>& inArray) {
    // Based on https://github.com/pybind/pybind11/issues/1377.
    // Request a buffer descriptor from Python.
    const auto buffer_info = inArray.request();

    // Extract data and shape of input array.
    T *data = static_cast<T *>(buffer_info.ptr);
    const auto shape = buffer_info.shape;

    // Wrap ndarray in Eigen::Map, which is then converted to an Eigen::Tensor.
    const Eigen::Tensor<T, 4>  in_tensor = Eigen::TensorMap<Eigen::Tensor<T, 4>>(data, shape[0], shape[1], shape[2], shape[3]);
    return in_tensor;
}


void bindCCSDEnvironment(py::module& module) {
    py::class_<GQCP::CCSDEnvironment<double>>(module, "CCSDEnvironment", "An algorithmic environment suitable for coupled-cluster calculations up to the CCSD level.")

        // CONSTRUCTORS
        .def_static(
            "PerturbativeCCSD",
            [](const GQCP::SQHamiltonian<double>& sq_hamiltonian, const GQCP::OrbitalSpace& orbital_space) {
                return GQCP::CCSDEnvironment<double>::PerturbativeCCSD(sq_hamiltonian, orbital_space);
            },
            "Initialize a CCSD algorithmic environment with initial guesses for the T1- and T2-amplitudes based on perturbation theory.")

        .def_static(
            "PerturbativeCCD",
            [](const GQCP::SQHamiltonian<double>& sq_hamiltonian, const GQCP::OrbitalSpace& orbital_space) {
                return GQCP::CCSDEnvironment<double>::PerturbativeCCD(sq_hamiltonian, orbital_space);
            },
            "Initialize a CCD algorithmic environment with initial guesses for the T2-amplitudes based on perturbation theory.")

        // Bind read-write members/properties, exposing intermediary environment variables to the Python interface.
        .def_readwrite("electronic_energies", &GQCP::CCSDEnvironment<double>::electronic_energies)


        // Define read-only 'getters'.
        .def_readonly(
            "t1_amplitudes",
            &GQCP::CCSDEnvironment<double>::t1_amplitudes)

        .def_readonly(
            "t2_amplitudes",
            &GQCP::CCSDEnvironment<double>::t2_amplitudes)


        // Bind methods for the replacement of the most current iterates.
        .def("replace_current_t1_amplitudes",
            [](GQCP::CCSDEnvironment<double>& environment, const GQCP::T1Amplitudes<double>& new_t1_amplitudes) {
                 environment.t1_amplitudes.pop_back();
                 environment.t1_amplitudes.push_back(new_t1_amplitudes);
             })

        .def("replace_current_t2_amplitudes",
            [](GQCP::CCSDEnvironment<double>& environment, const GQCP::T2Amplitudes<double>& new_t2_amplitudes) {
                 environment.t2_amplitudes.pop_back();
                 environment.t2_amplitudes.push_back(new_t2_amplitudes);
             })

        .def("replace_current_t2_amplitudes",
            [](GQCP::CCSDEnvironment<double>& environment, py::array_t<double>& new_t2_amplitudes, const GQCP::OrbitalSpace& orbital_space) {
                // Prepare the necessary members for ImplicitRankFourTensor.
                const auto axis1_indices = orbital_space.indices(GQCP::OccupationType::k_occupied);
                const auto axis2_indices = orbital_space.indices(GQCP::OccupationType::k_occupied);
                const auto axis3_indices = orbital_space.indices(GQCP::OccupationType::k_virtual);
                const auto axis4_indices = orbital_space.indices(GQCP::OccupationType::k_virtual);
                
                const auto t2_tensor = asTensor(new_t2_amplitudes);
                GQCP::T2Amplitudes<double> t2{GQCP::ImplicitRankFourTensorSlice<double>::FromIndices(axis1_indices, axis2_indices, axis3_indices, axis4_indices, t2_tensor), orbital_space};
                environment.t2_amplitudes.pop_back();
                environment.t2_amplitudes.push_back(t2);
            });
}


}  // namespace gqcpy

