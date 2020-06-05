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

#include "Mathematical/Grid/CubicGrid.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace py = pybind11;


namespace gqcpy {


void bindCubicGrid(py::module& module) {
    py::class_<GQCP::CubicGrid>(module, "CubicGrid", "A grid type whose points are on a regular cubic lattice.")

        // CONSTRUCTORS
        .def(py::init<>([](const Eigen::Vector3d& origin, const std::array<size_t, 3>& steps, const std::array<double, 3>& step_sizes) {
                 return GQCP::CubicGrid(GQCP::Vector<double, 3>(origin), steps, step_sizes);
             }),
             py::arg("origin"),
             py::arg("steps"),
             py::arg("step_sizes"))

        .def_static(
            "ReadCubeFile",
            &GQCP::CubicGrid::ReadCubeFile,
            py::arg("filename"),
            "Parse a GAUSSIAN Cube file (http://paulbourke.net/dataformats/cube/). The values for the contained scalar field are ignored.")

        .def_static(
            "ReadRegularGridFile",
            &GQCP::CubicGrid::ReadRegularGridFile,
            py::arg("filename"),
            "Parse an .rgrid-file and create the CubicGrid that is contained in it. The values for the scalar field or vector field are ignored.")


        // PUBLIC METHODS
        .def(
            "integrate",
            &GQCP::CubicGrid::integrate<double>,
            py::arg("field"),
            "Integrate a Field over this grid.")

        .def(
            "forEach",
            [](const GQCP::CubicGrid& cubic_grid, std::function<void(const size_t, const size_t, const size_t)>& callback) {
                cubic_grid.forEach(callback);
            },
            py::arg("callback"),
            "Loop over the points of this grid by index number.")

        .def(
            "forEach",
            [](const GQCP::CubicGrid& cubic_grid, const std::function<void(const Eigen::Vector3d&)>& callback) {
                const auto gqcp_callback = [callback](const GQCP::Vector<double, 3>& position) {
                    callback(Eigen::Vector3d(position));
                };

                cubic_grid.forEach(gqcp_callback);
            },
            py::arg("callback"),
            "Loop over the points of this grid by position (relative to the origin of this grid).")

        .def(
            "numberOfPoints",
            &GQCP::CubicGrid::numberOfPoints,
            "Return the number of points that are in this grid.")

        .def(
            "origin",
            [](const GQCP::CubicGrid& cubic_grid) {
                return Eigen::Vector3d(cubic_grid.origin());
            },
            "Return the origin of this grid.")

        .def(
            "position",
            [](const GQCP::CubicGrid& cubic_grid, size_t i, const size_t j, const size_t k) {
                return Eigen::Vector3d(cubic_grid.position(i, j, k));
            },
            py::arg("i"),
            py::arg("j"),
            py::arg("k"),
            "Return the position vector associated to the given indices.")

        .def(
            "numbersOfSteps",
            [](const GQCP::CubicGrid& cubic_grid, const size_t axis) {
                return cubic_grid.numbersOfSteps(axis);
            },
            py::arg("axis"),
            "Return the number of steps that can be taken in the direction of the specified axis.")

        .def(
            "numbersOfSteps",
            [](const GQCP::CubicGrid& cubic_grid) {
                return cubic_grid.numbersOfSteps();
            },
            "Return the number of steps in the x, y, z-directions.")

        .def(
            "stepSize",
            &GQCP::CubicGrid::stepSize,
            py::arg("axis"),
            "Return the step size that is taken in the direction of the specified axis.")

        .def(
            "stepSizes",
            &GQCP::CubicGrid::stepSizes,
            "Return the step sizes in the x, y, z-directions")

        .def(
            "totalVolume",
            &GQCP::CubicGrid::totalVolume,
            "Return the total volume that is contained in this grid.")

        .def(
            "writeToCubeFile",
            &GQCP::CubicGrid::writeToCubeFile,
            py::arg("scalar_field"),
            py::arg("filename"),
            py::arg("molecule"),
            "Write a field's values to a GAUSSIAN Cube file (http://paulbourke.net/dataformats/cube/).")

        .def(
            "voxelVolume",
            &GQCP::CubicGrid::voxelVolume,
            "Return the volume of the voxels in this grid.");
}


}  // namespace gqcpy
