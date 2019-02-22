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
#ifndef GQCP_NEWTONMINIMIZER_HPP
#define GQCP_NEWTONMINIMIZER_HPP



#include "math/optimization/BaseMinimizer.hpp"
#include "typedefs.hpp"



namespace GQCP {


/**
 *  An implementation of the Newton minimization algorithm
 */
class NewtonMinimizer : public BaseMinimizer {
private:
    VectorFunction grad;
    MatrixFunction H;


public:
    // CONSTRUCTORS
    /**
     *  @param x0                               the initial guess
     *  @param grad                             the callable gradient function
     *  @param H                                the callable Hessian function
     *  @param convergence_threshold            the threshold used for establishing convergence
     *  @param maximum_number_of_iterations     the maximum number of iterations in the algorithm
     */
    NewtonMinimizer(const Eigen::VectorXd& x0, const VectorFunction& grad, const MatrixFunction& H, double convergence_threshold = 1.0e-08, size_t maximum_number_of_iterations = 128);


    // DESTRUCTOR
    ~NewtonMinimizer() override = default;


    // PUBLIC OVERRIDDEN FUNCTIONS
    /**
     *  Minimize the function f(x)
     *
     *  If successful, sets
     *      - is_solved to true
     *      - the found solution
     */
    void solve() override;
};


}  // namespace GQCP



#endif  // GQCP_NEWTONMINIMIZER_HPP
