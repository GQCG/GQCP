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

#include "ONVBasis/ONVPath.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace py = pybind11;


namespace gqcpy {

void bindONVPath(py::module& module) {
    py::class_<GQCP::ONVPath>(module, "ONVPath", "A path-like representation of an ONV.")

        // CONSTRUCTOR

        .def(py::init<const GQCP::SpinUnresolvedONVBasis&, const GQCP::SpinUnresolvedONV&>(),
             py::arg("onv_basis"),
             py::arg("onv"))


        // PUBLIC METHODS

        .def(
            "address",
            [](const GQCP::ONVPath& path) {
                return path.address();
            },
            "Return the address of the current path")
        
        .def(
            "annihilate",
            [](GQCP::ONVPath& path) {
                path.annihilate();
            },
            "Annihilate the diagonal arc that starts at the current state of the path.")

        .def(
            "annihilate",
            [](GQCP::ONVPath& path, const size_t q, const size_t n) {
                path.annihilate(q, n);
            },
            py::arg("q"),
            py::arg("n"),
            "Annihilate the diagonal arc that starts at coordinate (q, n).")

        .def(
            "create",
            [](GQCP::ONVPath& path, const size_t p, const size_t n) {
                path.create(p, n);
            },
            py::arg("p"),
            py::arg("n"),
            "Create the diagonal arc that starts at coordinate (p, n). ")

        .def(
            "leftTranslate",
            [](GQCP::ONVPath& path, const size_t p, const size_t n) {
                path.leftTranslate(p, n);
            },
            py::arg("p"),
            py::arg("q"),
            "Translate the diagonal arc that starts at the coordinate (p, n) to the left.")

        .def(
            "nextCreationIndex",
            [](GQCP::ONVPath& path) {
                return path.nextCreationIndex();
            },
            "Return The orbital index p that should be checked next for a possible creation. Since we're always constructing paths from the top-left to the bottom-right, this index will always be larger than the index q on which we previously annihilated. After creation, the path then corresponds to E_{pq} |onv>, with |onv> the initial ONV.")

        .def(
            "leftTranslateUntilVertical",
            [](GQCP::ONVPath& path) {
                path.leftTranslateUntilVertical();
            },
            "Close the open path by shifting diagonal arcs to the left. Stop when an unoccupied orbital (vertical arc) is found.")

        .def(
            "sign",
            [](GQCP::ONVPath& path) {
                return path.sign();
            },
            "return the total phase factor/sign associated to the original path's modification");
}

}  // namespace gqcpy