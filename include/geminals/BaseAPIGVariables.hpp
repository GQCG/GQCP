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
#ifndef BaseAPIGVariables_hpp
#define BaseAPIGVariables_hpp


#include <Eigen/Dense>


namespace GQCP {


/**
 *  A base class for representing APIG-like variables.
 *
 *  In mathematical derivations, these are characterized by X_i^p, with i an occupied index, and p is a general orbital index. This class then effectively implements operator(i, p).
 */
class BaseAPIGVariables {
protected:
    size_t N_P;  // the number of electron pairs (= the number of geminals)
    size_t K;  // the number of orbitals
    Eigen::VectorXd x;  // the variables stored in a row-major form


public:
    // CONSTRUCTORS
    /**
     *  @param x        the variables in a vector representation that is in row-major storage
     *
     *  @param N_P      the number of electron pairs (= the number of geminals)
     *  @param K        the number of spatial orbitals
     */
    BaseAPIGVariables(const Eigen::VectorXd& x, size_t N_P, size_t K);

    /**
     *  Default constructor setting everything to zero
     */
    BaseAPIGVariables();  // default constructor needed


    // OPERATORS
    /**
     *  @param mu       a vector index
     *
     *  @return the variable x_mu
     */
    double operator()(size_t mu) const { return this->x(mu); }

    /**
     *  @param i        the major (geminal, subscript, non-contiguous) index
     *  @param p        the minor (orbital, superscript, contiguous) index
     *
     *  @return the variable X_i^p
     */
    double operator()(size_t i, size_t p) const;


    // GETTERS
    size_t get_N_P() const { return this->N_P; }
    size_t get_K() const { return this->K; }


    // PUBLIC METHODS
    /**
     *  @return the variables in row-major vector form
     */
    const Eigen::VectorXd& asVector() const { return this->x; }

    /**
     *  @return the variables in matrix form
     */
    virtual Eigen::MatrixXd asMatrix() const = 0;

    /**
     *  @param vector_index     the vector index of the variable
     *
     *  @return the major (geminal, subscript, non-contiguous) index i in the matrix of the variables
     */
    virtual size_t matrixIndexMajor(size_t vector_index) const = 0;

    /**
     *  @param vector_index     the vector index of the variable
     *
     *  @return the minor (orbital, superscript, contiguous) index p in the matrix of the variables
     */
    virtual size_t matrixIndexMinor(size_t vector_index) const = 0;

    /**
     *  @param i        the major (geminal, subscript, non-contiguous) index
     *  @param p        the minor (orbital, superscript, contiguous) index
     *
     *  @return the vector index of the variable X_i^p
     */
    virtual size_t vectorIndex(size_t i, size_t p) const = 0;
};


}  // namespace GQCP


#endif /* BaseAPIGVariables_hpp */
