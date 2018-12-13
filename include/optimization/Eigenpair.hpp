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
#ifndef NUMOPT_EIGENPAIR_HPP
#define NUMOPT_EIGENPAIR_HPP



#include <Eigen/Dense>



namespace numopt {
namespace eigenproblem {


/**
 *  A container class to store an eigenpair, i.e. an eigenvector with a corresponding eigenvalue
 */
class Eigenpair {
private:
    double eigenvalue;
    Eigen::VectorXd eigenvector;


public:
    // CONSTRUCTORS
    /**
     *  A constructor that sets the eigenvalue to zero and the corresponding eigenvector to zeros
     *
     *  @param dimension        the dimension of the eigenvector
     */
    explicit Eigenpair(size_t dimension = 1);

    /**
     *  @param eigenvalue       the eigenvalue
     *  @param eigenvector      the eigenvector
     */
    Eigenpair(double eigenvalue, const Eigen::VectorXd& eigenvector);


    // GETTERS
    double get_eigenvalue() const { return this->eigenvalue; };
    Eigen::VectorXd get_eigenvector() const { return this->eigenvector; };


    // PUBLIC METHODS
    /**
     *  @param other            the other Eigenpair
     *  @param tolerance        a tolerance for comparison
     *
     *  @return if this Eigenpair is equal to the other: if the eigenvalues and eigenvectors are equal given the tolerance
     */
    bool isEqual(const numopt::eigenproblem::Eigenpair& other, double tolerance=1.0e-08) const;
};


}  // namespace eigenproblem
}  // namespace numopt



#endif  // NUMOPT_EIGENPAIR_HPP
