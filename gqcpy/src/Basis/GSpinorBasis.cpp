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
#include "Basis/SpinorBasis/GSpinorBasis.hpp"

#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Molecule/Molecule.hpp"
#include "Operator/FirstQuantized/Operator.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace py = pybind11;


namespace gqcpy {


void bindGSpinorBasis(py::module& module) {
    py::class_<GQCP::GSpinorBasis<double, GQCP::GTOShell>>(module, "GSpinorBasis", "A class that represents a real, (generalized) spinor basis with underlying GTO shells.")

        .def(py::init<const GQCP::Molecule& , const std::string&>(), 
            py::arg("molecule"),
            py::arg("basisset_name")
        )

        .def("lowdinOrthonormalize",
            [ ] (GQCP::GSpinorBasis<double, GQCP::GTOShell>& spinor_basis) {
                spinor_basis.lowdinOrthonormalize();
            },
            "Transform the spinor basis to the 'Löwdin basis', which is the orthonormal basis that we transform to with T = S^{-1/2}, where S is the current overlap matrix."
        )

        .def("numberOfSpinors",
            [ ] (GQCP::GSpinorBasis<double, GQCP::GTOShell>& spinor_basis) {
                return spinor_basis.numberOfSpinors();
            },
            "Return the number of spinors that 'are' in this generalized spinor basis."
        )

        .def("quantizeKineticOperator",
            [ ] (const GQCP::GSpinorBasis<double, GQCP::GTOShell>& spinor_basis) {
                return spinor_basis.quantize(GQCP::Operator::Kinetic());
            },
            "Return the kinetic energy operator expressed in this spinor basis"
        )

        .def("quantizeNuclearAttractionOperator",
            [ ] (const GQCP::GSpinorBasis<double, GQCP::GTOShell>& spinor_basis, const GQCP::Molecule& molecule) {
                return spinor_basis.quantize(GQCP::Operator::NuclearAttraction(molecule));
            },
            "Return the nuclear attraction operator expressed in this spinor basis"
        )

        .def("quantizeOverlapOperator",
            [ ] (const GQCP::GSpinorBasis<double, GQCP::GTOShell>& spinor_basis) {
                return spinor_basis.quantize(GQCP::Operator::Overlap());
            },
            "Return the overlap operator expressed in this spinor basis"
        )

        .def("transform", [] (GQCP::GSpinorBasis<double, GQCP::GTOShell>& spinor_basis, const Eigen::MatrixXd& T_matrix) {
                const GQCP::TransformationMatrix<double> T (T_matrix);
                spinor_basis.transform(T);
            },
            "Transform the current spinor basis using a given transformation matrix"
        );
}


}  // namespace gqcpy
