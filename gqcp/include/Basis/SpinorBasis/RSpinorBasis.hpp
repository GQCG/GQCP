// This file is part of GQCG-GQCP.
//
// Copyright (C) 2017-2020  the GQCG developers
//
// GQCG-GQCP is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GQCG-GQCP is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-GQCP.  If not, see <http://www.gnu.org/licenses/>.

#pragma once


#include "Basis/Integrals/IntegralCalculator.hpp"
#include "Basis/ScalarBasis/ScalarBasis.hpp"
#include "Basis/SpinorBasis/GSpinorBasis.hpp"
#include "Basis/SpinorBasis/JacobiRotationParameters.hpp"
#include "Basis/SpinorBasis/SimpleSpinorBasis.hpp"
#include "Mathematical/Representation/QCMatrix.hpp"
#include "Molecule/Molecule.hpp"
#include "Molecule/NuclearFramework.hpp"
#include "Operator/FirstQuantized/Operator.hpp"
#include "Operator/SecondQuantized/SQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/SQTwoElectronOperator.hpp"
#include "Utilities/type_traits.hpp"
#include "Utilities/typedefs.hpp"

#include <Eigen/Dense>


namespace GQCP {


/**
 *  A class that represents a spinor basis in which the expansion of the alpha and beta components in terms of the underlying scalar orbitals are restricted to be equal
 * 
 *  @tparam _ExpansionScalar        the scalar type of the expansion coefficients
 *  @tparam _Shell                  the type of shell that the underlying scalar basis contains
 */
template <typename _ExpansionScalar, typename _Shell>
class RSpinorBasis:
    public SimpleSpinorBasis<_ExpansionScalar, RSpinorBasis<_ExpansionScalar, _Shell>> {

public:
    using Shell = _Shell;
    using BasisFunction = typename Shell::BasisFunction;
    using ExpansionScalar = _ExpansionScalar;

    using Base = SimpleSpinorBasis<_ExpansionScalar, RSpinorBasis<_ExpansionScalar, _Shell>>;


private:
    ScalarBasis<Shell> scalar_basis;  // the underlying scalar basis that is equal for both the alpha- and beta components


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param scalar_basis         the underlying scalar basis that is equal for both the alpha- and beta components
     *  @param C                    the matrix that holds the the expansion coefficients, i.e. that expresses the restricted spinors in terms of the underlying scalar basis
     */
    RSpinorBasis(const ScalarBasis<Shell>& scalar_basis, const TransformationMatrix<ExpansionScalar>& C) :
        Base(C),
        scalar_basis {scalar_basis} {}


    /**
     *  Construct a restricted spinor basis with an initial coefficient matrix that is the identity
     * 
     *  @param scalar_basis             the underlying scalar basis that is equal for both the alpha- and beta components
     * 
     *  @note the resulting restricted spinor basis is (most likely) non-orthogonal
     */
    RSpinorBasis(const ScalarBasis<Shell>& scalar_basis) :
        RSpinorBasis(scalar_basis, TransformationMatrix<double>::Identity(scalar_basis.numberOfBasisFunctions(),
                                                                          scalar_basis.numberOfBasisFunctions())) {}


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
        RSpinorBasis(ScalarBasis<Shell>(nuclear_framework, basisset_name)) {}


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
        RSpinorBasis(molecule.nuclearFramework(), basisset_name) {}


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @param ao_list     indices of the AOs used for the Mulliken populations
     *
     *  @return the Mulliken operator for a set of given AO indices
     *
     *  @note this method is only available for real matrix representations
     */
    template <typename Z = ExpansionScalar>
    enable_if_t<std::is_same<Z, double>::value, ScalarSQOneElectronOperator<double>> calculateMullikenOperator(const std::vector<size_t>& ao_list) const {

        const auto K = this->numberOfSpatialOrbitals();
        if (ao_list.size() > K) {
            throw std::invalid_argument("RSpinorBasis::calculateMullikenOperator(std::vector<size_t>): Too many AOs are selected in the given ao_list");
        }

        const SquareMatrix<double> p_a = SquareMatrix<double>::PartitionMatrix(ao_list, K);                       // the partitioning matrix
        const auto S_AO = IntegralCalculator::calculateLibintIntegrals(Operator::Overlap(), this->scalar_basis);  // the overlap matrix expressed in the AO basis

        ScalarSQOneElectronOperator<double> mulliken_op {0.5 * (this->C.adjoint() * p_a * S_AO * this->C + this->C.adjoint() * S_AO * p_a * this->C)};
        return mulliken_op;
    }


    /**
     *  @return this restricted spinor basis as a general one
     */
    GSpinorBasis<ExpansionScalar, Shell> generalized() const {

        // Build up the 'general' coefficient matrix.
        const auto K = this->numberOfSpatialOrbitals();
        const auto M = this->numberOfSpinors();
        TransformationMatrix<ExpansionScalar> C_general = TransformationMatrix<ExpansionScalar>::Zero(M, M);

        C_general.topLeftCorner(K, K) = this->coefficientMatrix();
        C_general.bottomRightCorner(K, K) = this->coefficientMatrix();

        return GSpinorBasis<ExpansionScalar, Shell>(this->scalarBasis(), C_general);  // the alpha- and beta- scalar bases are equal
    }


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
     *  @param fq_one_op                            the first-quantized one-electron operator
     * 
     *  @tparam FQOneElectronOperator               the type of the first-quantized one-electron operator
     * 
     *  @return the second-quantized operator corresponding to the given first-quantized operator
     */
    template <typename FQOneElectronOperator>
    auto quantize(const FQOneElectronOperator& fq_one_op) const -> SQOneElectronOperator<product_t<typename FQOneElectronOperator::Scalar, ExpansionScalar>, FQOneElectronOperator::Components> {

        using ResultScalar = product_t<typename FQOneElectronOperator::Scalar, ExpansionScalar>;
        using ResultOperator = SQOneElectronOperator<ResultScalar, FQOneElectronOperator::Components>;

        const auto one_op_par = IntegralCalculator::calculateLibintIntegrals(fq_one_op, this->scalarBasis());
        ResultOperator op {one_op_par};  // op for 'operator'
        op.transform(this->coefficientMatrix());
        return op;
    }


    /**
     *  @param fq_op        the first-quantized Coulomb operator
     * 
     *  @return the second-quantized operator corresponding to the Coulomb operator
     */
    auto quantize(const CoulombRepulsionOperator& fq_op) const -> SQTwoElectronOperator<product_t<CoulombRepulsionOperator::Scalar, ExpansionScalar>, CoulombRepulsionOperator::Components> {

        using ResultScalar = product_t<CoulombRepulsionOperator::Scalar, ExpansionScalar>;
        using ResultOperator = SQTwoElectronOperator<ResultScalar, CoulombRepulsionOperator::Components>;

        const auto one_op_par = IntegralCalculator::calculateLibintIntegrals(fq_op, this->scalarBasis());
        ResultOperator op {one_op_par};  // op for 'operator'
        op.transform(this->coefficientMatrix());
        return op;
    }


    /**
     *  @return the underlying scalar basis, which is equal for the alpha and beta components
     */
    const ScalarBasis<Shell>& scalarBasis() const { return this->scalar_basis; }
};


}  // namespace GQCP
