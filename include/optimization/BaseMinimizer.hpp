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
#ifndef NUMOPT_BASEMINIMIZER_HPP
#define NUMOPT_BASEMINIMIZER_HPP



#include <Eigen/Dense>



namespace numopt {
namespace minimization {


/**
 *  A base class for the implementation of minimizers
 */
class BaseMinimizer {
protected:
    constexpr static size_t maximum_number_of_iterations = 128;
    const double convergence_threshold = 1.0e-08;

    double is_solved = false;

    const Eigen::VectorXd x0;  // initial guess to the problem
    Eigen::VectorXd x;  // current guess or final solution to the problem


public:
    // CONSTRUCTORS
    /**
     *  Constructor based on a given initial guess @param x0 and a @param convergence_threshold
     */
    BaseMinimizer(const Eigen::VectorXd& x0, double convergence_threshold);


    // DESTRUCTOR
    virtual ~BaseMinimizer() = default;


    // GETTERS
    Eigen::VectorXd get_solution() const;


    // PUBLIC PURE VIRTUAL METHODS
    /**
     *  Solve the problem associated to the numerical minimization method
     *
     *  If successful, it sets
     *      - @member is_solved to true
     *      - @member x to the found solution
     */
    virtual void solve() = 0;
};


}  // namespace minimization
}  // namespace numopt



#endif  // NUMOPT_BASEMINIMIZER_HPP
