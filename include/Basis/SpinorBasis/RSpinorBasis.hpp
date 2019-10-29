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
 *  @tparam _ExpansionScalar        the scalar type of the expansion coefficients
 *  @tparam _Shell                  the type of shell that the underlying scalar basis contains
 */
template <typename _ExpansionScalar, typename _Shell>
class RSpinorBasis {
public:
    using Shell = _Shell;
    using BasisFunction = typename Shell::BasisFunction;
    using ExpansionScalar = _ExpansionScalar;


private:
    ScalarBasis<Shell> scalar_basis;  // the underlying scalar basis that is equal for both the alpha- and beta components
    TransformationMatrix<ExpansionScalar> C;  // the matrix that holds the the expansion coefficients, i.e. that expresses the restricted spinors in terms of the underlying scalar basis


public:

    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param scalar_basis         the underlying scalar basis that is equal for both the alpha- and beta components
     *  @param C                    the matrix that holds the the expansion coefficients, i.e. that expresses the restricted spinors in terms of the underlying scalar basis
     */
    RSpinorBasis(const ScalarBasis<Shell>& scalar_basis, const TransformationMatrix<ExpansionScalar>& C) :
        scalar_basis (scalar_basis),
        C (C)
    {}


    /**
     *  Construct a restricted spinor basis with an initial coefficient matrix that is the identity
     * 
     *  @param scalar_basis             the underlying scalar basis that is equal for both the alpha- and beta components
     * 
     *  @note the resulting restricted spinor basis is (most likely) non-orthogonal
     */
    RSpinorBasis(const ScalarBasis<Shell>& scalar_basis) : 
        RSpinorBasis(scalar_basis, TransformationMatrix<double>::Identity(scalar_basis.numberOfBasisFunctions(), scalar_basis.numberOfBasisFunctions()))
    {}


    /**
     *  Construct a restricted spinor basis with an underlying scalar basis (for both the alpha and beta components) by placing shells corresponding to the basisset specification on every nucleus of the nuclear framework
     *
     *  @param nuclear_framework        the nuclear framework containing the nuclei on which the shells should be centered
     *  @param basisset_name            the name of the basisset, e.g. "STO-3G"
     *
     *  @note the normalization factors of the spherical (or axis-aligned Cartesian) GTO primitives are embedded in the contraction coefficients of the underlying shells
     *  @note the resulting restricted spinor basis is (most likely) non-orthogonal
     */
    RSpinorBasis(const NuclearFramework& nuclear_framework, const std::string& basisset_name) :
        RSpinorBasis(ScalarBasis<Shell>(nuclear_framework, basisset_name))
    {}


    /**
     *  Construct a restricted spinor basis with an underlying scalar basis (for both the alpha and beta components) by placing shells corresponding to the basisset specification on every nucleus of the molecule
     *
     *  @param molecule             the molecule containing the nuclei on which the shells should be centered
     *  @param basisset_name        the name of the basisset, e.g. "STO-3G"
     *
     *  @note the normalization factors of the spherical (or axis-aligned Cartesian) GTO primitives are embedded in the contraction coefficients of the underlying shells
     *  @note the resulting restricted spinor basis is (most likely) non-orthogonal
     */
    RSpinorBasis(const Molecule& molecule, const std::string& basisset_name) :
        RSpinorBasis(molecule.nuclearFramework(), basisset_name)
    {}


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return the number of different spatial orbitals that are used in this restricted spinor basis
     */
    size_t numberOfSpatialOrbitals() const { return this->scalar_basis.numberOfBasisFunctions(); }

    /**
     *  @return the number of spinors that 'are' in this restricted spinor basis
     */
    size_t numberOfSpinors() const {
        return 2 * this->numberOfSpatialOrbitals();  // alpha and beta spinors are equal 
    }

    /**
     *  Transform the restricted spinor basis another one using the given transformation matrix
     * 
     *  @param T            the transformation matrix that transforms both the alpha- and beta components
     */
    void transform(const TransformationMatrix<ExpansionScalar>& T) {

        this->C.transform(T);
    }


    /**
     *  Rotate the restricted spinor basis to another one using the given unitary transformation matrix
     * 
     *  @param U            the unitary transformation matrix that transforms both the alpha- and beta components
     */
    void rotate(const TransformationMatrix<ExpansionScalar>& U) {

        // Check if the given matrix is actually unitary
        if (!U.isUnitary(1.0e-12)) {
            throw std::invalid_argument("RSpinorBasis::rotate(const TransformationMatrix<ExpansionScalar>&): The given transformation matrix is not unitary.");
        }

        this->transform(U);
    }


    /**
     *  Rotate the restricted spinor basis to another one using the unitary transformation matrix that corresponds to the given Jacobi rotation parameters
     * 
     *  @param jacobi_rotation_parameters       the Jacobi rotation parameters (p, q, angle) that are used to specify a Jacobi rotation: we use the (cos, sin, -sin, cos) definition for the Jacobi rotation matrix
     * 
     *  @note this function is only available for real restricted spinor bases because Jacobi rotation parameters generate real rotations
     */
    template<typename Z = ExpansionScalar>
    enable_if_t<std::is_same<Z, double>::value> rotate(const JacobiRotationParameters& jacobi_rotation_parameters) {

        const auto dim = this->numberOfSpatialOrbitals();
        const auto J = TransformationMatrix<double>::FromJacobi(jacobi_rotation_parameters, dim);
        this->rotate(J);
    }


    /**
     *  @return the transformation matrix between the scalar basis and the current orbitals
     */
    const TransformationMatrix<ExpansionScalar>& coefficientMatrix() const { return this->C; }

    /**
     *  @param precision                the precision used to test orthonormality
     * 
     *  @return if this restricted spinor basis basis is orthonormal within the given precision
     */
    bool isOrthonormal(const double precision = 1.0e-08) const {

        const auto S = this->overlapMatrix();

        const auto dim = this->numberOfSpatialOrbitals();
        return S.isApprox(SquareMatrix<ExpansionScalar>::Identity(dim, dim), precision);
    }


    /**
     *  Transform the restricted spinor basis to the Löwdin basis, which is the orthonormal basis that we transform to with T = S^{-1/2}, where S is the overlap matrix in the underlying scalar orbital basis
     */
    void lowdinOrthonormalize() {

        this->C = this->lowdinOrthonormalizationMatrix();
    }


    /**
     *  @return the transformation matrix to the Löwdin basis: T = S_current^{-1/2}
     */
    TransformationMatrix<double> lowdinOrthonormalizationMatrix() const {

        // Calculate S^{-1/2}, where S is epxressed in the current spinor basis
        const auto S = this->overlapMatrix();
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes (S);
        return TransformationMatrix<double>(saes.operatorInverseSqrt());
    }


    /**
     *  @return the overlap (one-electron) operator of this restricted spinor basis
     */
    ScalarSQOneElectronOperator<ExpansionScalar> overlap() const {

        return ScalarSQOneElectronOperator<ExpansionScalar>({this->overlapMatrix()});
    }


    /**
     *  @return the current overlap matrix of this restricted spinor basis
     */
    QCMatrix<ExpansionScalar> overlapMatrix() const {

        auto S = this->scalar_basis.calculateLibintIntegrals(Operator::Overlap());  // the overlap operator expressed in the underlying scalar basis
        S.basisTransformInPlace(this->C);  // the overlap operator expressed in this restricted spinor basis
        return S;
    }


    /**
     *  @return the underlying scalar basis, which is equal for the alpha and beta components
     */
    const ScalarBasis<Shell>& scalarBasis() const { return this->scalar_basis; }

    /**
     *  @param fq_op        the first-quantized operator
     * 
     *  @return the second-quantized operator corresponding to the given first-quantized operator
     */
    auto quantize(const OverlapOperator& fq_op) const -> SQOneElectronOperator<product_t<OverlapOperator::Scalar, ExpansionScalar>, OverlapOperator::Components> {

        using ResultScalar = product_t<OverlapOperator::Scalar, ExpansionScalar>;
        using ResultOperator = SQOneElectronOperator<ResultScalar, OverlapOperator::Components>;

        ResultOperator op ({this->scalarBasis().calculateLibintIntegrals(fq_op)});  // op for 'operator'
        op.transform(this->coefficientMatrix());
        return op;
    }


    /**
     *  @param fq_op        the first-quantized operator
     * 
     *  @return the second-quantized operator corresponding to the given first-quantized operator
     */
    auto quantize(const KineticOperator& fq_op) const -> SQOneElectronOperator<product_t<KineticOperator::Scalar, ExpansionScalar>, KineticOperator::Components> {

        using ResultScalar = product_t<KineticOperator::Scalar, ExpansionScalar>;
        using ResultOperator = SQOneElectronOperator<ResultScalar, KineticOperator::Components>;

        ResultOperator op ({this->scalarBasis().calculateLibintIntegrals(fq_op)});  // op for 'operator'
        op.transform(this->coefficientMatrix());
        return op;
    }


    /**
     *  @param fq_op        the first-quantized operator
     * 
     *  @return the second-quantized operator corresponding to the given first-quantized operator
     */
    auto quantize(const NuclearAttractionOperator& fq_op) const -> SQOneElectronOperator<product_t<NuclearAttractionOperator::Scalar, ExpansionScalar>, NuclearAttractionOperator::Components> {

        using ResultScalar = product_t<NuclearAttractionOperator::Scalar, ExpansionScalar>;
        using ResultOperator = SQOneElectronOperator<ResultScalar, NuclearAttractionOperator::Components>;

        ResultOperator op ({this->scalarBasis().calculateLibintIntegrals(fq_op)});  // op for 'operator'
        op.transform(this->coefficientMatrix());
        return op;
    }


    /**
     *  @param fq_op        the first-quantized operator
     * 
     *  @return the second-quantized operator corresponding to the given first-quantized operator
     */
    auto quantize(const ElectronicDipoleOperator& fq_op) const -> SQOneElectronOperator<product_t<ElectronicDipoleOperator::Scalar, ExpansionScalar>, ElectronicDipoleOperator::Components> {

        using ResultScalar = product_t<ElectronicDipoleOperator::Scalar, ExpansionScalar>;
        using ResultOperator = SQOneElectronOperator<ResultScalar, ElectronicDipoleOperator::Components>;

        ResultOperator op ({this->scalarBasis().calculateLibintIntegrals(fq_op)});  // op for 'operator'
        op.transform(this->coefficientMatrix());
        return op;
    }


    /**
     *  @param fq_op        the first-quantized operator
     * 
     *  @return the second-quantized operator corresponding to the given first-quantized operator
     */
    auto quantize(const CoulombRepulsionOperator& fq_op) const -> SQTwoElectronOperator<product_t<CoulombRepulsionOperator::Scalar, ExpansionScalar>, CoulombRepulsionOperator::Components> {

        using ResultScalar = product_t<CoulombRepulsionOperator::Scalar, ExpansionScalar>;
        using ResultOperator = SQTwoElectronOperator<ResultScalar, CoulombRepulsionOperator::Components>;

        ResultOperator op ({this->scalarBasis().calculateLibintIntegrals(fq_op)});  // op for 'operator'
        op.transform(this->coefficientMatrix());
        return op;
    }


    /**
     *  @param ao_list     indices of the AOs used for the Mulliken populations
     *
     *  @return the Mulliken operator for a set of AOs
     *
     *  @note this method is only available for real matrix representations
     */
    template<typename Z = ExpansionScalar>
    enable_if_t<std::is_same<Z, double>::value, ScalarSQOneElectronOperator<double>> calculateMullikenOperator(const Vectoru& ao_list) const {

        const auto K = this->numberOfSpatialOrbitals();
        if (ao_list.size() > K) {
            throw std::invalid_argument("RSpinorBasis::calculateMullikenOperator(Vectoru): Too many AOs are selected in the given ao_list");
        }

        const SquareMatrix<double> p_a = SquareMatrix<double>::PartitionMatrix(ao_list, K);  // the partitioning matrix
        const auto S_AO = this->scalar_basis.calculateLibintIntegrals(Operator::Overlap());  // the overlap matrix expressed in the AO basis

        ScalarSQOneElectronOperator<double> mulliken_op ({ (this->C.adjoint() * p_a * S_AO * this->C + this->C.adjoint() * S_AO * p_a * this->C)/2 });
        return mulliken_op;
    }
};


}  // namespace GQCP
