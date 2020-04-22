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


#include "Mathematical/Representation/BlockMatrix.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCModel/CI/LinearExpansion.hpp"
#include "QCModel/Geminals/GeminalCoefficientsInterface.hpp"


namespace GQCP {


/**
 *  Geminal coefficients for an AP1roG wave function.
 */
class AP1roGGeminalCoefficients: public GeminalCoefficientsInterface {
private:
    size_t K;               // the number of spatial orbitals corresponding to these geminal coefficients
    size_t N_P;             // the number of electron pairs (i.e. the number of geminals) corresponding to these geminal coefficients
    BlockMatrix<double> G;  // the AP1roG geminal coefficients (not including the identity matrix on the left), as a block matrix, so it implements easy operator(i,a) calls


public:
    // CONSTRUCTORS

    /**
     *  @param G            the AP1roG geminal coefficients (not including the identity matrix on the left), as a block matrix
     */
    AP1roGGeminalCoefficients(const BlockMatrix<double>& G, const size_t N_P, const size_t K);

    /**
     *  @param G            the AP1roG geminal coefficients (not including the identity matrix on the left)
     */
    AP1roGGeminalCoefficients(const MatrixX<double>& G);

    /**
     *  Constructor that sets the geminal coefficients to zero
     *
     *  @param N_P      the number of electron pairs (= the number of geminals)
     *  @param K        the number of spatial orbitals
     */
    AP1roGGeminalCoefficients(const size_t N_P, const size_t K);


    // DESTRUCTOR
    ~AP1roGGeminalCoefficients() override;


    // OPERATORS

    /**
     *  @param i            the zero-based index of the geminal, i.e. subscript of the geminal coefficient: i is in [0, N_P[ with N_P the number of electron pairs
     *  @param a            the zero-based index of the occupied orbital, i.e. superscript of the geminal coefficient: a is in [N_P, K[ with K the number of spatial orbitals
     * 
     *  @return an element of the AP1roG geminal coefficient matrix G_i^a
     */
    double operator()(const size_t i, const size_t a) const;


    // NAMED CONSTRUCTORS

    /**
     *  @param g        the geminal coefficients in a vector representation that is in row-major storage
     *
     *  @param N_P      the number of electron pairs (= the number of geminals)
     *  @param K        the number of spatial orbitals
     */
    static AP1roGGeminalCoefficients FromColumnMajor(const VectorX<double>& g, const size_t N_P, const size_t K);

    /**
     *  @param g        the geminal coefficients in a vector representation that is in row-major storage
     *
     *  @param N_P      the number of electron pairs (= the number of geminals)
     *  @param K        the number of spatial orbitals
     */
    static AP1roGGeminalCoefficients FromRowMajor(const VectorX<double>& g, const size_t N_P, const size_t K);

    /**
     *  @param sq_hamiltonian       the Hamiltonian expressed in an orthonormal basis
     *  @param N_P                  the number of electron pairs (= the number of geminals)
     *
     *  @return the AP1roG geminal coefficients in the weak interaction limit
     */
    static AP1roGGeminalCoefficients WeakInteractionLimit(const SQHamiltonian<double>& sq_hamiltonian, const size_t N_P);


    // STATIC PUBLIC METHODS

    /**
     *  @param N_P      the number of electron pairs (= the number of geminals)
     *  @param K        the number of spatial orbitals
     *
     *  @return the number of 'free' geminal coefficients
     */
    static size_t numberOfGeminalCoefficients(const size_t N_P, const size_t K);


    // GETTERS

    /**
     *  @return the number of electron pairs that these geminal coefficients describe
     */
    size_t get_N_P() const { return this->numberOfElectronPairs(); }

    /**
     *  @return the number of spatial orbitals each geminal is expanded in
     */
    size_t get_K() const { return this->numberOfSpatialOrbitals(); }


    // PUBLIC METHODS

    /**
     *  @return the total geminal coefficient matrix, including the identity matrix block
     */
    MatrixX<double> asMatrix() const;

    /**
     *  @return the geminal coefficients as a column-major vector, excluding the identity block
     */
    VectorX<double> asVector() const;

    /**
     *  @return the number of geminal coefficients that this instance encapsulates.
     */
    size_t count() const { return AP1roGGeminalCoefficients::numberOfGeminalCoefficients(this->N_P, this->K); }

    /**
     *  @return the number of electron pairs that these geminal coefficients describe
     */
    size_t numberOfElectronPairs() const { return this->N_P; }

    /**
     *  @return the number of spatial orbitals each geminal is expanded in
     */
    size_t numberOfSpatialOrbitals() const { return this->K; }

    /**
     *  @param onv      the doubly-occupied (spin-resolved) ONV that is being projected on
     *
     *  @return the overlap of the AP1roG wave function with the given ONV, i.e. the projection of the APIG wave function onto that ONV
     */
    double overlap(const SpinUnresolvedONV& onv) const override;
};


}  // namespace GQCP
