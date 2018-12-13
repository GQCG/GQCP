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
#ifndef NUMOPT_NEWTONMINIMIZER_HPP
#define NUMOPT_NEWTONMINIMIZER_HPP



#include "BaseMinimizer.hpp"
#include "common.hpp"



namespace numopt {
namespace minimization {


/**
 *  An implementation of the Newton minimization algorithm
 */
class NewtonMinimizer : public BaseMinimizer {
private:
    const VectorFunction grad;
    const MatrixFunction H;


public:
    // CONSTRUCTORS
    /**
     *  @param x0                           the initial guess
     *  @param grad                         the callable gradient function
     *  @param H                            the callable Hessian function
     *  @param convergence_threshold        the threshold used for establishing convergence
     */
    NewtonMinimizer(const Eigen::VectorXd& x0, const VectorFunction& grad, const MatrixFunction& H,
                    double convergence_threshold = 1.0e-08);


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


}  // namespace minimization
}  // numopt



#endif  // NUMOPT_NEWTONMINIMIZER_HPP
