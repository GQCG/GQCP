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

#include "Mathematical/Grid/Field.hpp"

#include <pybind11/eigen.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <functional>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


void bindField(py::module& module) {
    py::class_<Field<double>>(module, "ScalarField", "A set of function values corresponding to points in space.")

        // CONSTRUCTORS
        .def(py::init<const std::vector<double>&>(),
             py::arg("values"))

        .def_static(
            "ReadCubeFile",
            [](const std::string& filename) {
                return Field<double>::ReadCubeFile(filename);
            },
            py::arg("filename"),
            "Parse a GAUSSIAN Cube file (http://paulbourke.net/dataformats/cube/) for scalar field values. The grid-associated information is discarded.")


        // OPERATORS
        .def(py::self += py::self)

        .def(py::self + py::self)

        .def(-py::self)

        .def(py::self -= py::self)

        .def(py::self - py::self)


        // PUBLIC METHODS
        .def(
            "map",
            [](Field<double>& field, const std::function<double(const double&)>& function) {
                field.map(function);
            },
            py::arg("function"),
            "Apply a given function on each of this field's values, in-place.")

        .def(
            "mapped",
            [](const Field<double>& field, const std::function<double(const double&)>& function) {
                field.mapped(function);
            },
            py::arg("function"),
            "Apply a given function on each of this field's values.")

        .def(
            "size",
            &Field<double>::size,
            "Return the size of this field, i.e. the number of field values.")

        .def(
            "value",
            [](const Field<double>& field, const size_t index) {
                return field.value(index);
            },
            py::arg("index"),
            "Return a read-only field value, corresponding to the given index.")

        .def(
            "value",
            [](Field<double>& field, const size_t index) {
                return field.value(index);
            },
            py::arg("index"),
            "Return a read-only field value, corresponding to the given index.")

        .def(
            "values",
            &Field<double>::values,
            "Return the evaluated function values, in the order of the grid's loop.");


    py::class_<Field<Vector<double, 3>>>(module, "VectorField", "A set of function values corresponding to points in space.")

        // CONSTRUCTORS
        .def(py::init<>([](const std::vector<Eigen::Vector3d>& values) {
                 // Transform the Eigen::Vector3d vectors into Vector<double, 3>.
                 std::vector<Vector<double, 3>> gqcp_values;
                 gqcp_values.reserve(values.size());

                 std::transform(values.begin(), values.end(),
                                std::back_inserter(gqcp_values),
                                [](const Eigen::Vector3d& v) { return Vector<double, 3>(v); });

                 return Field<Vector<double, 3>>(gqcp_values);
             }),
             py::arg("values"))


        // OPERATORS
        .def(py::self += py::self)

        .def(py::self + py::self)

        .def(-py::self)

        .def(py::self -= py::self)

        .def(py::self - py::self)


        // PUBLIC METHODS
        .def(
            "map",
            [](Field<Vector<double, 3>>& field, const std::function<Vector<double, 3>(const Vector<double, 3>&)>& function) {
                field.map(function);
            },
            py::arg("function"),
            "Apply a given function on each of this field's values, in-place.")

        .def(
            "mapped",
            [](const Field<Vector<double, 3>>& field, const std::function<Vector<double, 3>(const Vector<double, 3>&)>& function) {
                field.mapped(function);
            },
            py::arg("function"),
            "Apply a given function on each of this field's values.")

        .def(
            "size",
            &Field<Vector<double, 3>>::size,
            "Return the size of this field, i.e. the number of field values.")

        .def(
            "value",
            [](const Field<Vector<double, 3>>& field, const size_t index) {
                return Eigen::Vector3d(field.value(index));
            },
            py::arg("index"),
            "Return a read-only field value, corresponding to the given index.")

        .def(
            "value",
            [](Field<Vector<double, 3>>& field, const size_t index) {
                return Eigen::Vector3d(field.value(index));
            },
            py::arg("index"),
            "Return a read-only field value, corresponding to the given index.")

        .def(
            "values",
            [](const Field<Vector<double, 3>>& field) {
                // Transform std::vector<Vector<double, 3>> into std::vector<Eigen::Vector3d>.
                const auto gqcp_values = field.values();

                std::vector<Eigen::Vector3d> eigen_values;
                eigen_values.reserve(gqcp_values.size());

                std::transform(gqcp_values.begin(), gqcp_values.end(),
                               std::back_inserter(eigen_values),
                               [](const Vector<double, 3>& v) {
                                   return Eigen::Vector3d(v);
                               });

                return eigen_values;
            },
            "Return the evaluated function values, in the order of the grid's loop.");
}


}  // namespace gqcpy
