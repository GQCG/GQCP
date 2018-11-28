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
#ifndef AP1roGGeminalCoefficients_hpp
#define AP1roGGeminalCoefficients_hpp


#include <Eigen/Dense>

#include "HamiltonianParameters/HamiltonianParameters.hpp"


namespace GQCP {


/**
 *  A class that represents geminal coefficients for AP1roG: the number of rows is the number of geminals, the number of columns is the number of spatial orbitals
 */
class AP1roGGeminalCoefficients {
private:
    size_t N_P;  // the number of electron pairs (= the number of geminals)
    size_t K;  // the number of orbitals
    Eigen::VectorXd g;  // the geminal coefficients stored in a row-major form


public:
    // CONSTRUCTORS
    /**
     *  Default constructor setting everything to zero
     */
    AP1roGGeminalCoefficients();
    
    /**
     *  Constructor that sets the geminal coefficients to zero
     *
     *  @param N_P      the number of electron pairs (= the number of geminals)
     *  @param K        the number of spatial orbitals
     */
    AP1roGGeminalCoefficients(size_t N_P, size_t K);

    /**
     *  @param g        the geminal coefficients in a vector representation that is in row-major storage
     *
     *  @param N_P      the number of electron pairs (= the number of geminals)
     *  @param K        the number of spatial orbitals
     */
    AP1roGGeminalCoefficients(const Eigen::VectorXd& g, size_t N_P, size_t K);

    /**
     *  @param ham_par      the Hamiltonian parameters
     *  @param N_P      the number of electron pairs (= the number of geminals)
     *
     *  @return the AP1roG geminal coefficients in the weak interaction limit
     */
    static AP1roGGeminalCoefficients WeakInteractionLimit(const HamiltonianParameters& ham_par, size_t N_P);


    // OPERATORS
    /**
     *  @param mu       a vector index
     *
     *  @return the geminal coefficient g_mu
     */
    double operator()(size_t mu) const;

    /**
     *  @param i        the major index (changes in i are not contiguous)
     *  @param a        the minor index (changes in a are contiguous)
     *
     *  @return the geminal coefficient G_i^a
     */
    double operator()(size_t i, size_t a) const;


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
     *  @param K        the number of spatial orbitals
     *  @param N_P      the number of electron pairs (= the number of geminals)
     *  @param vector_index     the vector index of the geminal coefficient
     *
     *  @return the major (non-contiguous) index i (i.e. the subscript) in the matrix of the geminal coefficients. Note that i is in [0 ... N_P[
     */
    static size_t matrixIndexMajor(size_t K, size_t N_P, size_t vector_index);

    /**
     *  @param K        the number of spatial orbitals
     *  @param N_P      the number of electron pairs (= the number of geminals)
     *  @param vector_index     the vector index of the geminal coefficient
     *
     *  @return the minor (contiguous) index a (i.e. the subscript) in the matrix of the geminal coefficients. Note that a is in [N_P ... K[
     */
    static size_t matrixIndexMinor(size_t K, size_t N_P, size_t vector_index);

    /**
     *  @param K        the number of spatial orbitals
     *  @param N_P      the number of electron pairs (= the number of geminals)
     *
     *  @param i        the major index (changes in i are not contiguous)
     *  @param a        the minor index (changes in a are contiguous)
     *
     *  @return the vector index of the geminal coefficient G_i^a
     */
    static size_t vectorIndex(size_t K, size_t N_P, size_t i, size_t a);


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
     *  @return the major (non-contiguous) index i (i.e. the subscript) in the matrix of the geminal coefficients. Note that i is in [0 ... N_P[
     */
    size_t matrixIndexMajor(size_t vector_index) const;

    /**
     *  @param vector_index     the vector index of the geminal coefficient
     *
     *  @return the minor (contiguous) index a (i.e. the subscript) in the matrix of the geminal coefficients. Note that a is in [N_P ... K[
     */
    size_t matrixIndexMinor(size_t vector_index) const;

    /**
     *  @param i        the major index (changes in i are not contiguous)
     *  @param a        the minor index (changes in a are contiguous)
     *
     *  @return the vector index of the geminal coefficient G_i^a
     */
    size_t vectorIndex(size_t i, size_t a) const;
};


}  // namespace GQCP



#endif /* AP1roGGeminalCoefficients_hpp */
