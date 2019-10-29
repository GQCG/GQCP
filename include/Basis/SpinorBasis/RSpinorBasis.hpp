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


#include "Basis/ScalarBasis/ScalarBasis.hpp"
#include "Basis/TransformationMatrix.hpp"
#include "Mathematical/Representation/QCMatrix.hpp"
#include "Molecule/Molecule.hpp"
#include "Molecule/NuclearFramework.hpp"
#include "Operator/FirstQuantized/Operator.hpp"
#include "Operator/SecondQuantized/SQOneElectronOperator.hpp"
#include "OrbitalOptimization/JacobiRotationParameters.hpp"

#include <Eigen/Dense>


namespace GQCP {


/**
 *  A class that represents a spinor basis in which the expansion of the alpha and beta components in terms of the underlying scalar orbitals are restricted to be equal
 * 
 *  @tparam _ShellType                  the type of shell that this scalar basis contains
 *  @tparam _TransformationScalar       the scalar type of the transformation matrix that connects the scalar basis with the current single-particle 'orbitals'
 */
template <typename _TransformationScalar, typename _ShellType>
class RSpinorBasis {
public:
    using ShellType = _ShellType;
    using BasisFunction = typename ShellType::BasisFunction;
    using TransformationScalar = _TransformationScalar;


private:
    ScalarBasis<ShellType> scalar_basis;  // the underlying scalar basis
    TransformationMatrix<TransformationScalar> T_total;  // the transformation matrix between the scalar basis and the current orbitals


public:

    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param scalar_basis             the underlying scalar basis
     *  @param T_total                  the transformation matrix between the scalar basis and the current orbitals
     */
    RSpinorBasis(const ScalarBasis<ShellType>& scalar_basis, const TransformationMatrix<TransformationScalar>& T_total) :
        scalar_basis (scalar_basis),
        T_total (T_total)
    {}


    /**
     *  Construct a single-particle basis with an identity transformation matrix
     * 
     *  @param scalar_basis             the underlying scalar basis
     */
    RSpinorBasis(const ScalarBasis<ShellType>& scalar_basis) : 
        RSpinorBasis(scalar_basis, TransformationMatrix<double>::Identity(scalar_basis.numberOfBasisFunctions(), scalar_basis.numberOfBasisFunctions()))
    {}


    /**
     *  Construct a single-particle basis by placing shells corresponding to the basisset specification on every nucleus of the nuclear framework
     *
     *  @param nuclear_framework        the nuclear framework containing the nuclei on which the shells should be centered
     *  @param basisset_name            the name of the basisset, e.g. "STO-3G"
     *
     *  @note the normalization factors of the spherical (or axis-aligned Cartesian) GTO primitives are embedded in the contraction coefficients of the underlying shells
     *  @note the resulting single-particle basis is (most likeley) non-orthogonal
     */
    RSpinorBasis(const NuclearFramework& nuclear_framework, const std::string& basisset_name) :
        RSpinorBasis(ScalarBasis<ShellType>(nuclear_framework, basisset_name))
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
    RSpinorBasis(const Molecule& molecule, const std::string& basisset_name) :
        RSpinorBasis(ScalarBasis<ShellType>(molecule, basisset_name))
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

        std::runtime_error("RSpinorBasis::basisFunctions(): This method hasn't been implemented yet");
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
    void transform(const TransformationMatrix<TransformationScalar>& T) {

        this->T_total.transform(T);
    }


    /**
     *  Rotate the single-particle basis to another one using the given unitary transformation matrix
     * 
     *  @param U            the unitary transformation matrix
     */
    void rotate(const TransformationMatrix<TransformationScalar>& U) {

        // Check if the given matrix is actually unitary
        if (!U.isUnitary(1.0e-12)) {
            throw std::invalid_argument("RSpinorBasis::rotate(const TransformationMatrix<TransformationScalar>&): The given transformation matrix is not unitary.");
        }

        this->transform(U);
    }


    /**
     *  Rotate the single-particle basis to another one using the unitary transformation matrix that corresponds to the given Jacobi rotation parameters
     * 
     *  @param jacobi_rotation_parameters       the Jacobi rotation parameters (p, q, angle) that are used to specify a Jacobi rotation: we use the (cos, sin, -sin, cos) definition for the Jacobi rotation matrix
     */
    void rotate(const JacobiRotationParameters& jacobi_rotation_parameters) {

        const auto dim = this->numberOfOrbitals();
        const auto J = TransformationMatrix<double>::FromJacobi(jacobi_rotation_parameters, dim);
        this->rotate(J);
    }


    /**
     *  @return the transformation matrix between the scalar basis and the current orbitals
     */
    TransformationMatrix<TransformationScalar> transformationMatrix() const { return this->T_total; }

    /**
     *  @param precision                the precision used to test orthonormality
     * 
     *  @return if this single-particle basis is orthonormal within the given precision
     */
    bool isOrthonormal(const double precision = 1.0e-08) const {

        const auto S = this->overlapMatrix();

        const auto dim = this->numberOfOrbitals();
        return S.isApprox(SquareMatrix<TransformationScalar>::Identity(dim, dim), precision);
    }


    /**
     *  Transform the single-particle basis to the Löwdin basis, which is the orthonormal basis that we transform to with T = S^{-1/2}, where S is the overlap matrix in the underlying scalar orbital basis
     */
    void lowdinOrthonormalize() {

        this->T_total = this->lowdinOrthonormalizationMatrix();
    }


    /**
     *  @return the transformation matrix to the Löwdin basis: T = S_current^{-1/2}
     */
    TransformationMatrix<double> lowdinOrthonormalizationMatrix() const {

        // The transformation matrix to the Löwdin basis is T = S_current^{-1/2}
        auto S = this->scalar_basis.calculateLibintIntegrals(Operator::Overlap());  // in the underlying (possibly orthonormal) scalar basis
        S.basisTransformInPlace(this->transformationMatrix());  // now S is expressed in the current orbital basis
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes (S);

        return TransformationMatrix<double>(saes.operatorInverseSqrt());
    }


    /**
     *  @return the overlap (one-electron) operator of this single-particle basis
     */
    ScalarSQOneElectronOperator<TransformationScalar> overlap() const {

        return ScalarSQOneElectronOperator<TransformationScalar>({this->overlapMatrix()});
    }



    /**
     *  @return the current overlap matrix of this single-particle basis
     */
    QCMatrix<TransformationScalar> overlapMatrix() const {

        auto S = this->scalar_basis.calculateLibintIntegrals(Operator::Overlap());  // in the underlying scalar basis
        S.basisTransformInPlace(this->T_total);  // in this single-particle basis
        return QCMatrix<TransformationScalar>(S);
    }



    /**
     *  @return the underlying scalar basis
     */
    const ScalarBasis<ShellType>& scalarBasis() const { return this->scalar_basis; }

    /**
     *  @param fq_op        the first-quantized operator
     * 
     *  @return the second-quantized operator corresponding to the given first-quantized operator
     */
    auto quantize(const OverlapOperator& fq_op) const -> SQOneElectronOperator<product_t<OverlapOperator::Scalar, TransformationScalar>, OverlapOperator::Components> {

        using ResultScalar = product_t<OverlapOperator::Scalar, TransformationScalar>;
        using ResultOperator = SQOneElectronOperator<ResultScalar, OverlapOperator::Components>;

        ResultOperator op ({this->scalarBasis().calculateLibintIntegrals(fq_op)});  // op for 'operator'
        op.transform(this->transformationMatrix());
        return op;
    }


    /**
     *  @param fq_op        the first-quantized operator
     * 
     *  @return the second-quantized operator corresponding to the given first-quantized operator
     */
    auto quantize(const KineticOperator& fq_op) const -> SQOneElectronOperator<product_t<KineticOperator::Scalar, TransformationScalar>, KineticOperator::Components> {

        using ResultScalar = product_t<KineticOperator::Scalar, TransformationScalar>;
        using ResultOperator = SQOneElectronOperator<ResultScalar, KineticOperator::Components>;

        ResultOperator op ({this->scalarBasis().calculateLibintIntegrals(fq_op)});  // op for 'operator'
        op.transform(this->transformationMatrix());
        return op;
    }


    /**
     *  @param fq_op        the first-quantized operator
     * 
     *  @return the second-quantized operator corresponding to the given first-quantized operator
     */
    auto quantize(const NuclearAttractionOperator& fq_op) const -> SQOneElectronOperator<product_t<NuclearAttractionOperator::Scalar, TransformationScalar>, NuclearAttractionOperator::Components> {

        using ResultScalar = product_t<NuclearAttractionOperator::Scalar, TransformationScalar>;
        using ResultOperator = SQOneElectronOperator<ResultScalar, NuclearAttractionOperator::Components>;

        ResultOperator op ({this->scalarBasis().calculateLibintIntegrals(fq_op)});  // op for 'operator'
        op.transform(this->transformationMatrix());
        return op;
    }


    /**
     *  @param fq_op        the first-quantized operator
     * 
     *  @return the second-quantized operator corresponding to the given first-quantized operator
     */
    auto quantize(const ElectronicDipoleOperator& fq_op) const -> SQOneElectronOperator<product_t<ElectronicDipoleOperator::Scalar, TransformationScalar>, ElectronicDipoleOperator::Components> {

        using ResultScalar = product_t<ElectronicDipoleOperator::Scalar, TransformationScalar>;
        using ResultOperator = SQOneElectronOperator<ResultScalar, ElectronicDipoleOperator::Components>;

        ResultOperator op ({this->scalarBasis().calculateLibintIntegrals(fq_op)});  // op for 'operator'
        op.transform(this->transformationMatrix());
        return op;
    }


    /**
     *  @param fq_op        the first-quantized operator
     * 
     *  @return the second-quantized operator corresponding to the given first-quantized operator
     */
    auto quantize(const CoulombRepulsionOperator& fq_op) const -> SQTwoElectronOperator<product_t<CoulombRepulsionOperator::Scalar, TransformationScalar>, CoulombRepulsionOperator::Components> {

        using ResultScalar = product_t<CoulombRepulsionOperator::Scalar, TransformationScalar>;
        using ResultOperator = SQTwoElectronOperator<ResultScalar, CoulombRepulsionOperator::Components>;

        ResultOperator op ({this->scalarBasis().calculateLibintIntegrals(fq_op)});  // op for 'operator'
        op.transform(this->transformationMatrix());
        return op;
    }


    /**
     *  @param ao_list     indices of the AOs used for the Mulliken populations
     *
     *  @return the Mulliken operator for a set of AOs
     *
     *  Note that this method is only available for real matrix representations
     */
    template<typename Z = TransformationScalar>
    enable_if_t<std::is_same<Z, double>::value, ScalarSQOneElectronOperator<double>> calculateMullikenOperator(const Vectoru& ao_list) const {

        const auto K = this->numberOfBasisFunctions();
        if (ao_list.size() > K) {
            throw std::invalid_argument("RSpinorBasis::calculateMullikenOperator(Vectoru): Too many AOs are selected");
        }

        // Create the partitioning matrix
        SquareMatrix<double> p_a = SquareMatrix<double>::PartitionMatrix(ao_list, K);

        ScalarSQOneElectronOperator<double> S_AO ({this->overlapMatrix()});
        TransformationMatrix<double> T_inverse = this->transformationMatrix().inverse();
        S_AO.transform(T_inverse);

        ScalarSQOneElectronOperator<double> mulliken_matrix ({ (T_total.adjoint() * p_a * S_AO.parameters() * T_total + T_total.adjoint() * S_AO.parameters() * p_a * T_total)/2 });

        return mulliken_matrix;
    }
};


}  // namespace GQCP
