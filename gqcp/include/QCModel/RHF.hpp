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


#include "Basis/TransformationMatrix.hpp"
#include "Processing/RDM/OneRDM.hpp"
#include "QCMethod/QCObjective.hpp"


namespace GQCP {
namespace QCModel {


/**
 *  The restricted Hartree-Fock wave function model
 * 
 *  @tparam _ExpansionScalar_       the type of scalar that is used for the expansion of the spatial orbitals in their underlying scalar basis
 */
template <typename _ExpansionScalar>
class RHF {
public:
    using ExpansionScalar = _ExpansionScalar;


private:
    size_t N_P;  // the number of electron pairs
    TransformationMatrix<ExpansionScalar> C;  // the coefficient matrix that expresses every spatial orbital (as a column) in its underlying scalar basis


public:
    // STATIC PUBLIC METHODS

    /**
     *  @param C    the coefficient matrix that expresses every spatial orbital (as a column) in its underlying scalar basis
     *  @param N    the number of electrons
     *
     *  @return the RHF 1-RDM expressed in the underlying scalar basis
     */
    static OneRDM<double> calculateScalarBasis1RDM(const TransformationMatrix<double>& C, const size_t N) {

        const size_t K = C.dimension();
        const auto D_orthonormal = RHF::calculateOrthonormalBasis1RDM(K, N);

        // Transform the 1-RDM in an orthonormal basis to the underlying scalar basis
        return C.conjugate() * D_orthonormal * C.transpose();
    }

    /**
     *  @param K    the number of spatial orbitals
     *  @param N    the number of electrons
     *
     *  @return the RHF 1-RDM expressed in an orthonormal spinor basis
     */
    static OneRDM<ExpansionScalar> calculateOrthonormalBasis1RDM(const size_t K, const size_t N) {

        if (N % 2 != 0) {
            throw std::invalid_argument("RHF::calculateOrthonormalBasis1RDM(const size_t, const size_t): The number of given electrons cannot be odd for RHF.");
        }

        // The 1-RDM for RHF looks like (for K=5, N=6)
        //    2  0  0  0  0
        //    0  2  0  0  0
        //    0  0  2  0  0
        //    0  0  0  0  0
        //    0  0  0  0  0

        OneRDM<double> D_MO = OneRDM<double>::Zero(K, K);
        D_MO.topLeftCorner(N/2, N/2) = 2 * SquareMatrix<double>::Identity(N/2, N/2);

        return D_MO;
    }

    // PUBLIC METHODS

    /**
     *  @return the coefficient matrix that expresses every spatial orbital (as a column) in its underlying scalar basis
     */
    const TransformationMatrix<ExpansionScalar>& coefficientMatrix() const { return this->C; }

    /**
     *  @return an objective that can check if the RHF Fock matrix is diagonal
     */
    static QCObjective<RHF> DiagonalFockMatrix();

    /**
     *  @return the number of spatial orbitals that these RHF model parameters describe
     */
    size_t numberOfSpatialOrbitals() const { this->coefficientMatrix.dimension(); }
};


}  // namespace QCModel
}  // namespace GQCP
