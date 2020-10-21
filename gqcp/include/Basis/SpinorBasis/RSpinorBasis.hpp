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
#include "Basis/SpinorBasis/SimpleSpinorBasis.hpp"
#include "Basis/SpinorBasis/Spinor.hpp"
#include "Basis/Transformations/JacobiRotationParameters.hpp"
#include "Basis/Transformations/RTransformationMatrix.hpp"
#include "Mathematical/Representation/SquareMatrix.hpp"
#include "Molecule/Molecule.hpp"
#include "Molecule/NuclearFramework.hpp"
#include "Operator/FirstQuantized/Operator.hpp"
#include "Operator/SecondQuantized/EvaluatableScalarRSQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/SQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/SQTwoElectronOperator.hpp"
#include "Utilities/aliases.hpp"
#include "Utilities/type_traits.hpp"

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
    // The scalar type used to represent an expansion coefficient of the spinors in the underlying scalar orbitals: real or complex.
    using ExpansionScalar = _ExpansionScalar;

    // The type of shell the underlying scalar bases contain.
    using Shell = _Shell;

    // The type of the base spinor basis.
    using Base = SimpleSpinorBasis<_ExpansionScalar, RSpinorBasis<_ExpansionScalar, _Shell>>;

    // The type of transformation matrix that is naturally related to a GSpinorBasis.
    using TM = RTransformationMatrix<ExpansionScalar>;  // TODO: Rename to TransformationMatrix once the class is gone

    using Primitive = typename Shell::Primitive;
    using BasisFunction = typename Shell::BasisFunction;
    using SpatialOrbital = LinearCombination<product_t<ExpansionScalar, typename BasisFunction::CoefficientScalar>, BasisFunction>;


    // The type of one-electron operators that are naturally associated with scalar operators expressed in this restricted spin-orbital bases.
    // using ScalarSQOneElectronOperator = ScalarRSQOneElectronOperator<ExpansionScalar>;

    // The type of one-electron operators that are naturally associated with scalar operators expressed in this restricted spin-orbital bases.
    // using ScalarSQTwoElectronOperator = ScalarRSQTwoElectronOperator<ExpansionScalar>;


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
    RSpinorBasis(const ScalarBasis<Shell>& scalar_basis, const TM& C) :
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
        RSpinorBasis(scalar_basis, TM::Identity(scalar_basis.numberOfBasisFunctions())) {}


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
    template <typename S = ExpansionScalar, typename = IsReal<S>>
    ScalarRSQOneElectronOperator<double> calculateMullikenOperator(const std::vector<size_t>& ao_list) const {

        const auto K = this->numberOfSpatialOrbitals();
        if (ao_list.size() > K) {
            throw std::invalid_argument("RSpinorBasis::calculateMullikenOperator(std::vector<size_t>): Too many AOs are selected in the given ao_list");
        }

        const SquareMatrix<double> p_a = SquareMatrix<double>::PartitionMatrix(ao_list, K);                       // the partitioning matrix
        const auto S_AO = IntegralCalculator::calculateLibintIntegrals(Operator::Overlap(), this->scalar_basis);  // the overlap matrix expressed in the AO basis

        ScalarRSQOneElectronOperator<double> mulliken_op {0.5 * (this->C.adjoint() * p_a * S_AO * this->C + this->C.adjoint() * S_AO * p_a * this->C)};
        return mulliken_op;
    }

    /**
     *  @return the number of different spatial orbitals that are used in this restricted spinor basis
     */
    size_t numberOfSpatialOrbitals() const { return this->scalar_basis.numberOfBasisFunctions(); }

    /**
     *  @return the number of spinors that 'are' in this restricted spinor basis
     */
    size_t numberOfSpinors() const { return 2 * this->numberOfSpatialOrbitals(); /* alpha and beta spinors are equal*/ }


    /**
     *  @param fq_op        the first-quantized Coulomb operator
     * 
     *  @return the second-quantized operator corresponding to the Coulomb operator
     */
    auto quantize(const CoulombRepulsionOperator& fq_op) const -> SQTwoElectronOperator<product_t<CoulombRepulsionOperator::Scalar, ExpansionScalar>, CoulombRepulsionOperator::NumberOfComponents> {

        using ResultScalar = product_t<CoulombRepulsionOperator::Scalar, ExpansionScalar>;
        using ResultOperator = SQTwoElectronOperator<ResultScalar, CoulombRepulsionOperator::NumberOfComponents>;

        const auto one_op_par = IntegralCalculator::calculateLibintIntegrals(fq_op, this->scalarBasis());  // in AO/scalar basis

        ResultOperator op {one_op_par};           // op for 'operator'
        op.transform(this->coefficientMatrix());  // now in spatial/spin-orbital basis
        return op;
    }


    /**
     *  Quantize the (one-electron) electronic density operator.
     * 
     *  @param fq_density_op                    the first-quantized density operator
     * 
     *  @return the second-quantized density operator
     */
    EvaluatableScalarRSQOneElectronOperator<ScalarFunctionProduct<LinearCombination<double, BasisFunction>>> quantize(const ElectronicDensityOperator& fq_density_op) const {

        using Evaluatable = ScalarFunctionProduct<LinearCombination<double, BasisFunction>>;  // the evaluatable type for the density operator
        using ResultOperator = EvaluatableScalarRSQOneElectronOperator<Evaluatable>;

        // There aren't any 'integrals' to be calculated for the density operator: we can just multiply every pair of spatial orbitals.
        const auto phi = this->spatialOrbitals();
        const auto K = this->numberOfSpatialOrbitals();

        SquareMatrix<Evaluatable> rho_par {K};  // the matrix representation ('par' for 'parameters') of the second-quantized (one-electron) density operator
        for (size_t p = 0; p < K; p++) {
            for (size_t q = 0; q < K; q++) {
                rho_par(p, q) = phi[p] * phi[q];
            }
        }

        return ResultOperator(rho_par);
    }


    /**
     *  @param fq_one_op                            the first-quantized one-electron operator
     * 
     *  @tparam FQOneElectronOperator               the type of the first-quantized one-electron operator
     * 
     *  @return the second-quantized operator corresponding to the given first-quantized operator
     */
    template <typename FQOneElectronOperator>
    auto quantize(const FQOneElectronOperator& fq_one_op) const -> RSQOneElectronOperator<product_t<typename FQOneElectronOperator::Scalar, ExpansionScalar>, typename FQOneElectronOperator::Vectorizer> {

        using ResultScalar = product_t<typename FQOneElectronOperator::Scalar, ExpansionScalar>;
        using ResultOperator = RSQOneElectronOperator<ResultScalar, typename FQOneElectronOperator::Vectorizer>;

        const auto one_op_par = IntegralCalculator::calculateLibintIntegrals(fq_one_op, this->scalarBasis());  // in AO/scalar basis

        ResultOperator op {one_op_par};           // op for 'operator'
        op.transform(this->coefficientMatrix());  // now in the spatial/spin-orbital basis
        return op;
    }


    /**
     *  @return the underlying scalar basis, which is equal for the alpha and beta components
     */
    const ScalarBasis<Shell>& scalarBasis() const { return this->scalar_basis; }

    /**
     *  @return the set of spatial orbitals that is associated to this spin-orbital basis
     */
    std::vector<SpatialOrbital> spatialOrbitals() const {

        // The spatial orbitals are a linear combination of the basis functions, where every column of the coefficient matrix describes one expansion of a spatial orbital in terms of the basis functions.
        const auto basis_functions = this->scalar_basis.basisFunctions();
        const auto& C = this->C;


        // For all spatial orbitals, proceed to calculate the contraction between the associated coefficient matrix column and the basis functions.
        std::vector<SpatialOrbital> spatial_orbitals;
        spatial_orbitals.reserve(this->numberOfSpatialOrbitals());
        for (size_t p = 0; p < this->numberOfSpatialOrbitals(); p++) {

            // Calculate the spatial orbitals as a contraction between a column of the coefficient matrix and the basis functions.
            SpatialOrbital spatial_orbital {};
            for (size_t mu = 0; mu < basis_functions.size(); mu++) {
                const auto coefficient = this->C.col(p)(mu);
                const auto& function = basis_functions[mu];
                spatial_orbital.append({coefficient}, {function});
            }

            spatial_orbitals.push_back(spatial_orbital);
        }

        return spatial_orbitals;
    }


    /**
     *  @return the set of spin-orbitals that is associated to this spin-orbital basis
     */
    std::vector<Spinor<ExpansionScalar, BasisFunction>> spinOrbitals() const {

        // The spin-orbitals for a restricted spin-orbital basis can be easily constructed from the spatial orbitals, by assigning a zero component once for the beta component of the spin-orbital and once for the alpha component of the spin-orbital.
        const auto spatial_orbitals = this->spatialOrbitals();

        std::vector<Spinor<ExpansionScalar, BasisFunction>> spin_orbitals;
        spin_orbitals.reserve(this->numberOfSpinors());
        for (const auto& spatial_orbital : spatial_orbitals) {

            // Add the alpha- and beta-spin-orbitals accordingly.
            const Spinor<ExpansionScalar, BasisFunction> alpha_spin_orbital {spatial_orbital, 0};  // the '0' int literal can be converted to a zero LinearCombination
            spin_orbitals.push_back(alpha_spin_orbital);

            const Spinor<ExpansionScalar, BasisFunction> beta_spin_orbital {0, spatial_orbital};
            spin_orbitals.push_back(beta_spin_orbital);
        }

        return spin_orbitals;
    }
};


/*
 *  MARK: SpinorBasisTraits
 */

/**
 *  A type that provides compile-time information on spinor bases that is otherwise not accessible through a public class alias.
 */
template <typename _ExpansionScalar, typename _Shell>
class SpinorBasisTraits<RSpinorBasis<_ExpansionScalar, _Shell>> {
public:
    // The scalar type used to represent an expansion coefficient of the spinors in the underlying scalar orbitals: real or complex.
    using ExpansionScalar = _ExpansionScalar;

    // The type of transformation matrix that is naturally related to an RSpinorBasis.
    using TM = RTransformationMatrix<ExpansionScalar>;  // TODO: Rename to TransformationMatrix once the class is gone

    // The second-quantized representation of the overlap operator related to an RSpinorBasis.
    using SQOverlapOperator = ScalarRSQOneElectronOperator<ExpansionScalar>;
};


}  // namespace GQCP
