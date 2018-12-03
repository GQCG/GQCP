// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2018  the GQCG developers
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
#ifndef BaseAPIGGeminalCoefficients_hpp
#define BaseAPIGGeminalCoefficients_hpp


#include <Eigen/Dense>

#include "WaveFunction/WaveFunction.hpp"


namespace GQCP {


/**
 *  A base class for representing APIG-like geminal wave functions
 */
class BaseAPIGGeminalCoefficients {
protected:
    size_t N_P;  // the number of electron pairs (= the number of geminals)
    size_t K;  // the number of orbitals
    Eigen::VectorXd g;  // the geminal coefficients stored in a row-major form


public:
    // CONSTRUCTORS
    /**
     *  @param g        the geminal coefficients in a vector representation that is in row-major storage
     *
     *  @param N_P      the number of electron pairs (= the number of geminals)
     *  @param K        the number of spatial orbitals
     */
    BaseAPIGGeminalCoefficients(const Eigen::VectorXd& g, size_t N_P, size_t K);

    /**
     *  Default constructor setting everything to zero
     */
    BaseAPIGGeminalCoefficients();


    // OPERATORS
    /**
     *  @param mu       a vector index
     *
     *  @return the geminal coefficient g_mu
     */
    double operator()(size_t mu) const { return this->g(mu); }

    /**
     *  @param i        the major (geminal, subscript, non-contiguous) index
     *  @param p        the minor (orbital, superscript, contiguous) index
     *
     *  @return the geminal coefficient G_i^p
     */
    double operator()(size_t i, size_t p) const;


    // GETTERS
    size_t get_N_P() const { return this->N_P; }
    size_t get_K() const { return this->K; }


    // PUBLIC METHODS
    /**
     *  @return the geminal coefficients in row-major vector form
     */
    const Eigen::VectorXd& asVector() const { return this->g; }

    /**
     *  @return the geminal coefficients in matrix form
     */
    virtual Eigen::MatrixXd asMatrix() const = 0;

    /**
     *  @param vector_index     the vector index of the geminal coefficient
     *
     *  @return the major (geminal, subscript, non-contiguous) index i in the matrix of the geminal coefficients
     */
    virtual size_t matrixIndexMajor(size_t vector_index) const = 0;

    /**
     *  @param vector_index     the vector index of the geminal coefficient
     *
     *  @return the minor (orbital, superscript, contiguous) index p in the matrix of the geminal coefficients
     */
    virtual size_t matrixIndexMinor(size_t vector_index) const = 0;

    /**
     *  @param i        the major (geminal, subscript, non-contiguous) index
     *  @param p        the minor (orbital, superscript, contiguous) index
     *
     *  @return the vector index of the geminal coefficient G_i^p
     */
    virtual size_t vectorIndex(size_t i, size_t p) const = 0;

    /**
     *  @param onv      the ONV that is being projected on
     *
     *  @return the overlap of the APIG-like wave function with the given on, i.e. the projection of the APIG wave function onto that ONV
     */
    virtual double overlap(const ONV& onv) const = 0;

    /**
     *  @return the wave function expansion corresponding to the geminal coefficients
     */
    WaveFunction toWaveFunction() const;
};


}  // namespace GQCP


#endif /* BaseAPIGGeminalCoefficients_hpp */
