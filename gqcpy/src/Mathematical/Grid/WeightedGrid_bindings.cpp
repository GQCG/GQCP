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

#include "Mathematical/Grid/WeightedGrid.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace py = pybind11;


namespace gqcpy {


void bindWeightedGrid(py::module& module) {
    py::class_<GQCP::WeightedGrid>(module, "WeightedGrid", "A collection of points in 3D-space, with each point associated to a weight.")

        // CONSTRUCTORS
        .def(py::init<>([](const std::vector<Eigen::Vector3d>& points, const Eigen::ArrayXd& weights) {
                 // Transform the Eigen::Vector3d vectors into GQCP::Vector<double, 3>.
                 std::vector<GQCP::Vector<double, 3>> gqcp_points;
                 gqcp_points.reserve(points.size());

                 std::transform(points.begin(), points.end(),
                                std::back_inserter(gqcp_points),
                                [](const Eigen::Vector3d& v) { return GQCP::Vector<double, 3>(v); });

                 return GQCP::WeightedGrid(gqcp_points, GQCP::ArrayX<double>(weights));
             }),
             py::arg("points"),
             py::arg("weights"))

        .def_static(
            "ReadIntegrationGridFile",
            &GQCP::WeightedGrid::ReadIntegrationGridFile,
            py::arg("filename"),
            "Parse an .igrid-file and create the WeightedGrid that is contained in it. The values for the scalar field or vector field are discarded.")


        // PUBLIC METHODS
        .def(
            "integrate",
            &GQCP::WeightedGrid::integrate<double>,
            py::arg("field"),
            "Integrate a Field over this grid.")

        .def(
            "numberOfPoints",
            &GQCP::WeightedGrid::numberOfPoints,
            "Return the number of points that are in this grid.")

        .def(
            "point",
            [](const GQCP::WeightedGrid& weighted_grid, const size_t index) {
                return Eigen::Vector3d(weighted_grid.point(index));
            },
            py::arg("index"),
            "Access one of the grid's points.")

        .def(
            "points",
            [](const GQCP::WeightedGrid& weighted_grid) {
                // Transform std::vector<GQCP::Vector<double, 3>> into std::vector<Eigen::Vector3d>.
                const auto gqcp_points = weighted_grid.points();

                std::vector<Eigen::Vector3d> eigen_points;
                eigen_points.reserve(gqcp_points.size());

                std::transform(gqcp_points.begin(), gqcp_points.end(),
                               std::back_inserter(eigen_points),
                               [](const GQCP::Vector<double, 3>& v) {
                                   return Eigen::Vector3d(v);
                               });

                return eigen_points;
            },
            "Return the grid points.")

        .def(
            "size",
            &GQCP::WeightedGrid::size,
            "Return the size of the grid, i.e. the number of grid points/weights.")

        .def(
            "weight",
            &GQCP::WeightedGrid::weight,
            py::arg("index"),
            "Return a read-only grid weight, corresponding to the given index.")

        .def(
            "weights",
            [](const GQCP::WeightedGrid& weighted_grid) {
                return Eigen::ArrayXd(weighted_grid.weights());
            },
            "Return a 1-D array containing the weights for each of the grid points.");
}


}  // namespace gqcpy
