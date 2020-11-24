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

#include "Basis/Transformations/GTransformation.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace py = pybind11;


namespace gqcpy {

using namespace GQCP;


void bindGTransformation(py::module& module) {

    py::class_<GTransformation<double>>(module, "GTransformation", "A 'general' basis transformation, i.e. a general, full-spinor basis transformation where the transformation mixes the alpha- and beta components of the two-component spinors.")


        /*
         *  MARK: Constructors
         */

        .def(py::init<>([](const Eigen::MatrixXd& T) {
                 return GTransformation<double> {T};
             }),
             py::arg("T"))


        /*
         *  MARK: Transformation matrix
         */

        .def("matrix",
             &GTransformation<double>::matrix,
             "Return the transformation matrix that collects the expansion coefficients of the new basis (vectors) in the old basis as columns.");
}


}  // namespace gqcpy
