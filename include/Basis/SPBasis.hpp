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


#include "Basis/ScalarBasis.hpp"
#include "Mathematical/Representation/SquareMatrix.hpp"
#include "OrbitalOptimization/JacobiRotationParameters.hpp"

#include <Eigen/Dense>


namespace GQCP {


/**
 *  A class that represents a single-particle basis
 * 
 *  @tparam _ShellType                  the type of shell that this scalar basis contains
 *  @tparam _TransformationScalar       the scalar type of the transformation matrix that connects the scalar basis with the current single-particle 'orbitals'
 */
template <typename _TransformationScalar, typename _ShellType>
class SPBasis {
public:
    using ShellType = _ShellType;
    using BasisFunction = typename ShellType::BasisFunction;
    using TransformationScalar = _TransformationScalar;


private:
    ScalarBasis<ShellType> scalar_basis;  // the underlying scalar basis
    SquareMatrix<TransformationScalar> T_total;  // the transformation matrix between the scalar basis and the current orbitals


public:

    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param scalar_basis             the underlying scalar basis
     *  @param T_total                  the transformation matrix between the scalar basis and the current orbitals
     */
    SPBasis(const ScalarBasis<ShellType>& scalar_basis, const SquareMatrix<TransformationScalar>& T_total) :
        scalar_basis (scalar_basis),
        T_total (T_total)
    {}


    /**
     *  Construct a single-particle basis with an identity transformation matrix
     * 
     *  @param scalar_basis             the underlying scalar basis
     */
    SPBasis(const ScalarBasis<ShellType>& scalar_basis) : 
        SPBasis(scalar_basis, SquareMatrix<double>::Identity(scalar_basis.numberOfBasisFunctions(), scalar_basis.numberOfBasisFunctions()))
    {}

    /**
     *  Construct a single-particle basis by placing shells corresponding to the basisset specification on every nucleus of the molecule
     *
     *  @param molecule             the molecule containing the nuclei on which the shells should be centered
     *  @param basisset_name        the name of the basisset, e.g. "STO-3G"
     *
     *  @note the normalization factors of the spherical (or axis-aligned Cartesian) GTO primitives are embedded in the contraction coefficients of the underlying shells
     *  @note the resulting single-particle basis is (most likeley) non-orthogonal
     */
    SPBasis(const Molecule& molecule, const std::string& basisset_name) :
        SPBasis(ScalarBasis<ShellType>(molecule, basisset_name))
    {}


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return the number of orbitals that 'are' in this single-particle basis
     */
    size_t numberOfBasisFunctions() const { return this->scalar_basis.numberOfBasisFunctions(); }

    /**
     *  @return the number of orbitals that 'are' in this single-particle basis
     */
    size_t numberOfOrbitals() const { return this->numberOfBasisFunctions(); }

    /**
     *  @return the orbitals that 'are' in this single-particle basis
     */
    std::vector<LinearCombination<double, BasisFunction>> basisFunctions() const { 

        std::runtime_error("SPBasis::basisFunctions(): This method hasn't been implemented yet");
    }


    /**
     *  @param i            the index of the requested orbital
     * 
     *  @return the orbital with the given index that 'is' in this scalar basis
     */
    LinearCombination<double, BasisFunction> basisFunction(const size_t i) const { return this->basisFunctions()[i]; }

    /**
     *  Transform the single-particle basis to another one using the given transformation matrix
     * 
     *  @param T            the transformation matrix
     */
    void transform(const SquareMatrix<TransformationScalar>& T) {

        this->T_total = this->T_total * T;
    }


    /**
     *  Rotate the single-particle basis to another one using the given unitary transformation matrix
     * 
     *  @param U            the unitary transformation matrix
     */
    void rotate(const SquareMatrix<TransformationScalar>& U) {

        // Check if the given matrix is actually unitary
        if (!U.isUnitary(1.0e-12)) {
            throw std::invalid_argument("SPBasis::rotate(const SquareMatrix<TransformationScalar>&): The given transformation matrix is not unitary.");
        }

        this->transform(U);
    }


    /**
     *  Rotate the single-particle basis to another one using the unitary transformation matrix that corresponds to the given Jacobi rotation parameters
     * 
     *  @param jacobi_rotation_parameters       the Jacobi rotation parameters (p, q, angle) that are used to specify a Jacobi rotation: we use the (cos, sin, -sin, cos) definition for the Jacobi rotation matrix. See transform() for how the transformation matrix between the two bases should be represented
     */
    void rotate(const JacobiRotationParameters& jacobi_rotation_parameters) {

        const auto dim = this->numberOfOrbitals();
        const auto J = SquareMatrix<double>::FromJacobi(jacobi_rotation_parameters, dim);
        this->rotate(J);
    }


    /**
     *  @return the transformation matrix between the scalar basis and the current orbitals
     */
    SquareMatrix<TransformationScalar> transformationMatrix() const { return this->T_total; }

    /**
     *  @param precision                the precision used to test orthonormality
     * 
     *  @return if this single-particle basis is orthonormal within the given precision
     */
    bool isOrthonormal(const double precision = 1.0e-08) const {

        const auto S = this->overlapMatrix();

        const auto dim = this->numberOfOrbitals();
        return S.isApprox(SquareMatrix<TransformationScalar>::Identity(this->dim, this->dim), precision);
    }


    /**
     *  Transform the single-particle basis to the Löwdin basis, which is the orthonormal basis that we transform to with T = S^{-1/2}, where S is the overlap matrix in the underlying scalar orbital basis
     * 
     *  @return transformation matrix to the Löwdin basis: T = S^{-1/2}
     */
    SquareMatrix<TransformationScalar> LowdinOrthonormalize() {

        // The transformation matrix to the Löwdin basis is T = S^{-1/2}
        auto S = this->scalar_basis.calculateLibintOverlapIntegrals();  // in the underlying (possibly orthonormal) scalar basis
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes (S);

        this->T_total = SquareMatrix<double>(saes.operatorInverseSqrt());
    }


    /**
     *  @return the current overlap matrix of this single-particle basis
     */
    SquareMatrix<TransformationScalar> overlapMatrix() const {

        auto S = this->scalar_basis.calculateLibintOverlapIntegrals();  // in the underlying scalar basis
        S.basisTransformInPlace(this->T_total);  // in this single-particle basis
        return S;
    }
};


}  // namespace GQCP
