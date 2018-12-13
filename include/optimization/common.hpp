// This file is part of GQCG-numopt.
// 
// Copyright (C) 2017-2018  the GQCG developers
// 
// GQCG-numopt is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// GQCG-numopt is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-numopt.  If not, see <http://www.gnu.org/licenses/>.
#ifndef NUMOPT_COMMON_HPP
#define NUMOPT_COMMON_HPP

#include <Eigen/Dense>


/*
 *  TYPEDEFS
 */

namespace numopt {


using VectorFunction = std::function<Eigen::VectorXd (const Eigen::VectorXd&)>;
using MatrixFunction = std::function<Eigen::MatrixXd (const Eigen::VectorXd&)>;


}  // namespace numopt



#endif  // NUMOPT_COMMON_HPP
