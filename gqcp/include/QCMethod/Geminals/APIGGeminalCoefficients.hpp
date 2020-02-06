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
#pragma once


#include "ONVBasis/LinearExpansion/LinearExpansion.hpp"
#include "Mathematical/Representation/Matrix.hpp"
#include "QCMethod/Geminals/GeminalCoefficientsInterface.hpp"



namespace GQCP {


/**
 *  A class that represents geminal coefficients for an APIG wave function
 */
class APIGGeminalCoefficients : public GeminalCoefficientsInterface {
private:
    size_t K;  // the number of spatial orbitals corresponding to these geminal coefficients
    size_t N_P;  // the number of electron pairs (i.e. the number of geminals) corresponding to these geminal coefficients
    MatrixX<double> G;  // the APIG geminal coefficients


public:
    // CONSTRUCTORS

    /**
     *  @param G                the APIG geminal coefficients in their matrix representation
     */
    APIGGeminalCoefficients(const MatrixX<double>& G);

    /**
     *  Constructor that sets the geminal coefficients to zero
     *
     *  @param N_P      the number of electron pairs (= the number of geminals)
     *  @param K        the number of spatial orbitals
     */
    APIGGeminalCoefficients(const size_t N_P, const size_t K);

    // /**
    //  *  Default constructor setting everything to zero
    //  */
    // APIGGeminalCoefficients();  // default constructor needed


    // DESTRUCTOR
    ~APIGGeminalCoefficients() override;


    // OPERATORS

    /**
     *  @param i            the zero-based index of the geminal, i.e. subscript of the geminal coefficient: i is in [0, N_P[ with N_P the number of electron pairs
     *  @param p            the zero-based index of the orbital, i.e. superscript of the geminal coefficient
     * 
     *  @return an element of the AP1roG geminal coefficient matrix G_i^a
     */
    double operator()(const size_t i, const size_t p) const;


    // NAMED CONSTRUCTORS

    /**
     *  @param g        the geminal coefficients in a vector representation that is in column-major storage
     *  @param N_P      the number of electron pairs (= the number of geminals)
     *  @param K        the number of spatial orbitals
     */
    static APIGGeminalCoefficients FromColumnMajor(const VectorX<double>& g, const size_t N_P, const size_t K);

    /**
     *  @param g        the geminal coefficients in a vector representation that is in row-major storage
     *  @param N_P      the number of electron pairs (= the number of geminals)
     *  @param K        the number of spatial orbitals
     */
    static APIGGeminalCoefficients FromRowMajor(const VectorX<double>& g, const size_t N_P, const size_t K);


    // STATIC PUBLIC METHODS

    /**
     *  @param N_P      the number of electron pairs (= the number of geminals)
     *  @param K        the number of spatial orbitals
     *
     *  @return the number of 'free' geminal coefficients
     */
    static size_t numberOfGeminalCoefficients(size_t N_P, size_t K);


    // PUBLIC METHODS

    /**
     *  @return the geminal coefficients in matrix form
     */
    const MatrixX<double>& asMatrix() const { return this->G; };

    /**
     *  @param onv      the ONV that is being projected on
     *
     *  @return the overlap of the APIG wave function with the given on, i.e. the projection of the APIG wave function onto that ONV
     */
    double overlap(const ONV& onv) const override;
};


}  // namespace GQCP
