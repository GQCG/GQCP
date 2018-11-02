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
     *  @ return the geminal coefficient G_i^a
     */
    double operator()(size_t i, size_t a) const;


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
    Eigen::MatrixXd asMatrix() const;

    /**
     *  @param N_P      the number of electron pairs (= the number of geminals)
     *  @param K        the number of spatial orbitals
     *
     *  @return the number of 'free' geminal coefficients
     */
    static size_t numberOfGeminalCoefficients(size_t N_P, size_t K);

    /**
     *  For a geminal coefficient g_mu, return its major index in the matrix of geminal coefficients.
     *
     *      Note that:
     *          - the major index is i (i.e. the subscript), since changes in i are not contiguous
     *          - i is in [0 ... N_P[
     */
    size_t matrixIndexMajor(size_t vector_index) const;

    /**
     *  For a geminal coefficient g_mu, return its minor index in the matrix of geminal coefficients.
     *
     *      Note that:
     *          - the minor index is a (i.e. the superscript), since changes in a are contiguous
     *          - a is in [N_P ... K[
     */
    size_t matrixIndexMinor(size_t vector_index) const;

    /**
     *  For a geminal coefficient G_i^a, return its index in the vector of geminal coefficients.
     *
     *      Note that
     *          - i is in [0 ... N_P[       is the 'major' index (i.e. changes in i are not contiguous)
     *          - a is in [N_P ... K[       is the 'minor' index (i.e. changes in a are contiguous)
     */
    size_t vectorIndex(size_t i, size_t a) const;
};


}  // namespace GQCP



#endif /* AP1roGGeminalCoefficients_hpp */
