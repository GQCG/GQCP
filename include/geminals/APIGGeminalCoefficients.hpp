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
#ifndef APIGGeminalCoefficients_hpp
#define APIGGeminalCoefficients_hpp


#include <Eigen/Dense>

#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "WaveFunction/WaveFunction.hpp"



namespace GQCP {


/**
 *  A class that represents geminal coefficients for an APIG wave function
 */
class APIGGeminalCoefficients {
private:
    size_t N_P;  // the number of electron pairs (= the number of geminals)
    size_t K;  // the number of orbitals
    Eigen::VectorXd g;  // the geminal coefficients stored in a row-major form


public:
    // CONSTRUCTORS
    /**
     *  Constructor that sets the geminal coefficients to zero
     *
     *  @param N_P      the number of electron pairs (= the number of geminals)
     *  @param K        the number of spatial orbitals
     */
    APIGGeminalCoefficients(size_t N_P, size_t K);

    /**
     *  @param g        the geminal coefficients in a vector representation that is in row-major storage
     *
     *  @param N_P      the number of electron pairs (= the number of geminals)
     *  @param K        the number of spatial orbitals
     */
    APIGGeminalCoefficients(const Eigen::VectorXd& g, size_t N_P, size_t K);


    // OPERATORS
    /**
     *  @param mu       a vector index
     *
     *  @return the geminal coefficient g_mu
     */
    double operator()(size_t mu) const;

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


    // STATIC PUBLIC METHODS
    /**
     *  @param N_P      the number of electron pairs (= the number of geminals)
     *  @param K        the number of spatial orbitals
     *
     *  @return the number of 'free' geminal coefficients
     */
    static size_t numberOfGeminalCoefficients(size_t N_P, size_t K);

    /**
     *  @param K                the number of spatial orbitals
     *  @param N_P              the number of electron pairs (= the number of geminals)
     *  @param vector_index     the vector index of the geminal coefficient
     *
     *  @return the major (geminal, subscript, non-contiguous) index i in the matrix of the geminal coefficients
     */
    static size_t matrixIndexMajor(size_t K, size_t N_P, size_t vector_index);

    /**
     *  @param K                the number of spatial orbitals
     *  @param N_P              the number of electron pairs (= the number of geminals)
     *  @param vector_index     the vector index of the geminal coefficient
     *
     *  @return the minor (orbital, superscript, contiguous) index p in the matrix of the geminal coefficients
     */
    static size_t matrixIndexMinor(size_t K, size_t N_P, size_t vector_index);

    /**
     *  @param K        the number of spatial orbitals
     *  @param N_P      the number of electron pairs (= the number of geminals)
     *
     *  @param i        the major (geminal, subscript, non-contiguous) index
     *  @param p        the minor (orbital, superscript, contiguous) index
     *
     *  @return the vector index of the geminal coefficient G_i^p
     */
    static size_t vectorIndex(size_t K, size_t N_P, size_t i, size_t p);


    // PUBLIC METHODS
    /**
     *  @return the geminal coefficients in row-major vector form
     */
    const Eigen::VectorXd& asVector() const { return this->g; }

    /**
     *  @return the geminal coefficients in matrix form
     */
    Eigen::MatrixXd asMatrix() const;

    /**
     *  @param vector_index     the vector index of the geminal coefficient
     *
     *  @return the major (geminal, subscript, non-contiguous) index i in the matrix of the geminal coefficients
     */
    size_t matrixIndexMajor(size_t vector_index) const;

    /**
     *  @param vector_index     the vector index of the geminal coefficient
     *
     *  @return the minor (orbital, superscript, contiguous) index p in the matrix of the geminal coefficients
     */
    size_t matrixIndexMinor(size_t vector_index) const;

    /**
     *  @param i        the major (geminal, subscript, non-contiguous) index
     *  @param p        the minor (orbital, superscript, contiguous) index
     *
     *  @return the vector index of the geminal coefficient G_i^p
     */
    size_t vectorIndex(size_t i, size_t p) const;

    /**
     *  @return the wave function expansion corresponding to the geminal coefficients
     */
    WaveFunction toWaveFunction() const;
};


}  // namespace GQCP



#endif /* APIGGeminalCoefficients_hpp */
