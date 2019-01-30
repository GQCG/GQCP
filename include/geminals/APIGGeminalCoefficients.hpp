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
#ifndef APIGGeminalCoefficients_hpp
#define APIGGeminalCoefficients_hpp


#include <Eigen/Dense>

#include "geminals/BaseAPIGVariables.hpp"
#include "geminals/GeminalCoefficientsInterface.hpp"
#include "WaveFunction/WaveFunction.hpp"



namespace GQCP {


/**
 *  A class that represents geminal coefficients for an APIG wave function
 */
class APIGGeminalCoefficients : public BaseAPIGVariables, public GeminalCoefficientsInterface {
public:
    // CONSTRUCTORS
    /**
     *  Default constructor setting everything to zero
     */
    APIGGeminalCoefficients();  // default constructor needed

    /**
     *  @param g        the geminal coefficients in a vector representation that is in row-major storage
     *
     *  @param N_P      the number of electron pairs (= the number of geminals)
     *  @param K        the number of spatial orbitals
     */
    APIGGeminalCoefficients(const Eigen::VectorXd& g, size_t N_P, size_t K);

    /**
     *  Constructor that sets the geminal coefficients to zero
     *
     *  @param N_P      the number of electron pairs (= the number of geminals)
     *  @param K        the number of spatial orbitals
     */
    APIGGeminalCoefficients(size_t N_P, size_t K);

    /**
     *  @param G        the geminal coefficients in a matrix representation
     */
    APIGGeminalCoefficients(const Eigen::MatrixXd& G);


    // DESTRUCTOR
    ~APIGGeminalCoefficients() override;


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
    Eigen::MatrixXd asMatrix() const override;

    /**
     *  @param vector_index     the vector index of the geminal coefficient
     *
     *  @return the major (geminal, subscript, non-contiguous) index i in the matrix of the geminal coefficients
     */
    size_t matrixIndexMajor(size_t vector_index) const override;

    /**
     *  @param vector_index     the vector index of the geminal coefficient
     *
     *  @return the minor (orbital, superscript, contiguous) index p in the matrix of the geminal coefficients
     */
    size_t matrixIndexMinor(size_t vector_index) const override;

    /**
     *  @param i        the major (geminal, subscript, non-contiguous) index
     *  @param p        the minor (orbital, superscript, contiguous) index
     *
     *  @return the vector index of the geminal coefficient G_i^p
     */
    size_t vectorIndex(size_t i, size_t p) const override;

    /**
     *  @param onv      the ONV that is being projected on
     *
     *  @return the overlap of the APIG wave function with the given on, i.e. the projection of the APIG wave function onto that ONV
     */
    double overlap(const ONV& onv) const override;
};


}  // namespace GQCP



#endif /* APIGGeminalCoefficients_hpp */
