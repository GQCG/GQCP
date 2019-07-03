
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "Atom.hpp"
#include "Drivers/HubbardDriver.hpp"
#include "Drivers/FullConfigurationInteractionDriver.hpp"

namespace py = pybind11;

PYBIND11_MODULE(gqcpy, m) {

    py::class_<GQCP::HubbardDriver>(m, "HubbardDriver")
            .def(py::init<std::string, size_t, size_t, size_t, size_t>())
            .def("get_energies", &GQCP::HubbardDriver::get_energies)
            .def("get_first_order_rdms", &GQCP::HubbardDriver::get_first_order_rdms);

    py::class_<GQCP::FullConfigurationInteractionDriver>(m, "FullConfigurationInteractionDriver")
            .def(py::init<std::string, std::string, size_t, size_t>())
            .def("get_energy", &GQCP::FullConfigurationInteractionDriver::get_energy);
}
