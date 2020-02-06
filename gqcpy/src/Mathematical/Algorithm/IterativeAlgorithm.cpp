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
#include "Mathematical/Algorithm/IterativeAlgorithm.hpp"
#include "QCMethod/HF/RHFSCFEnvironment.hpp"

#include <pybind11/pybind11.h>


namespace py = pybind11;


namespace gqcpy {


void bindRHFSCFIterativeAlgorithm(py::module& module) {
    py::class_<GQCP::IterativeAlgorithm<GQCP::RHFSCFEnvironment<double>>>(module, "RHFSCFIterativeAlgorithm", "An algorithm that performs iterations using an RHFSCFEnvironment. In every iteration, convergence is checked and a set of iteration steps is performed.");
}


}  // namespace gqcpy
