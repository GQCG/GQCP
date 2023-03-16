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


#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Basis/ScalarBasis/LondonGTOShell.hpp"
#include "Basis/SpinorBasis/SimpleSpinOrbitalBasis.hpp"
#include "Basis/Transformations/JacobiRotation.hpp"
#include "Basis/Transformations/UTransformationComponent.hpp"
#include "Domain/UMullikenDomainComponent.hpp"
#include "Operator/SecondQuantized/USQOneElectronOperatorComponent.hpp"
#include "Utilities/type_traits.hpp"


namespace GQCP {


/**
 *  A type specifically designed as one of the spin-components of a `USpinOrbitalBasis`.
 * 
 *  @tparam _ExpansionScalar        The scalar type used to represent an expansion coefficient of the spin-orbitals in the underlying scalar orbitals: real or complex.
 *  @tparam _Shell                  The type of shell the underlying scalar basis contains.
 */
template <typename _ExpansionScalar, typename _Shell>
class USpinOrbitalBasisComponent:
    public SimpleSpinOrbitalBasis<_ExpansionScalar, _Shell, USpinOrbitalBasisComponent<_ExpansionScalar, _Shell>> {
public:
    // The scalar type used to represent an expansion coefficient of the spinors in the underlying scalar orbitals: real or complex.
    using ExpansionScalar = _ExpansionScalar;

    // The type of shell the underlying scalar bases contain.
    using Shell = _Shell;

    // The type that is used for representing the primitive for a basis function of this spin-orbital basis' underlying AO basis.
    using Primitive = typename Shell::Primitive;

    // The type that is used for representing the underlying basis functions of this spin-orbital basis.
    using BasisFunction = typename Shell::BasisFunction;


public:
    /*
     *  MARK: Constructors
     */

    // Inherit `SimpleSpinOrbitalBasis`'s constructors.
    using SimpleSpinOrbitalBasis<ExpansionScalar, Shell, USpinOrbitalBasisComponent<ExpansionScalar, Shell>>::SimpleSpinOrbitalBasis;


    /*
     *  MARK: Quantization of first-quantized operators (GTOShell)
     */

    /**
     *  Quantize a first-quantized one-electron operator.
     * 
     *  @param fq_one_op                            The first-quantized one-electron operator.
     * 
     *  @tparam FQOneElectronOperator               The type of the first-quantized one-electron operator.
     * 
     *  @return The second-quantized operator corresponding to the given first-quantized operator, i.e. expressed in/projected onto this spin-orbital basis.
     */
    template <typename FQOneElectronOperator, typename Z = Shell>
    auto quantize(const FQOneElectronOperator& fq_one_op) const -> enable_if_t<std::is_same<Z, GTOShell>::value, USQOneElectronOperatorComponent<product_t<typename FQOneElectronOperator::Scalar, ExpansionScalar>, typename FQOneElectronOperator::Vectorizer>> {

        // FIXME: Try to move this API to SimpleSpinOrbitalBasis.
        using ResultScalar = product_t<typename FQOneElectronOperator::Scalar, ExpansionScalar>;
        using ResultOperator = USQOneElectronOperatorComponent<ResultScalar, typename FQOneElectronOperator::Vectorizer>;

        const auto f_par = IntegralCalculator::calculateLibintIntegrals(fq_one_op, this->scalarBasis());  // Currently f is expressed in the AO basis. 'par' for 'parameters'

        ResultOperator f {f_par};
        f.transform(this->expansion());  // Now, f is expressed in the spin-orbital basis.
        return f;
    }


    /*
     *  MARK: Quantization of first-quantized operators (LondonGTOShell)
     */

    /**
     *  Quantize a first-quantized one-electron operator.
     * 
     *  @param fq_one_op                            The first-quantized one-electron operator.
     * 
     *  @tparam FQOneElectronOperator               The type of the first-quantized one-electron operator.
     * 
     *  @return The second-quantized operator corresponding to the given first-quantized operator, i.e. expressed in/projected onto this spin-orbital basis.
     */
    template <typename FQOneElectronOperator, typename Z = Shell>
    auto quantize(const FQOneElectronOperator& fq_one_op) const -> enable_if_t<std::is_same<Z, LondonGTOShell>::value, USQOneElectronOperatorComponent<product_t<typename FQOneElectronOperator::Scalar, ExpansionScalar>, typename FQOneElectronOperator::Vectorizer>> {

        using ResultScalar = product_t<typename FQOneElectronOperator::Scalar, ExpansionScalar>;
        using ResultOperator = USQOneElectronOperatorComponent<ResultScalar, typename FQOneElectronOperator::Vectorizer>;
        using Vectorizer = typename FQOneElectronOperator::Vectorizer;

        const auto N = FQOneElectronOperator::NumberOfComponents;
        const auto& vectorizer = FQOneElectronOperator::vectorizer;


        // The strategy for calculating the matrix representation of the one-electron operator in this spin-orbital basis is to:
        //  1. Express the operator in the underlying scalar bases and;
        //  2. Afterwards transform them using the current coefficient matrix.
        const auto K = this->scalarBasis().numberOfBasisFunctions();


        // 1. Express the operator in the underlying scalar basis.
        auto engine = GQCP::IntegralEngine::InHouse<GQCP::LondonGTOShell>(fq_one_op);
        const auto F = GQCP::IntegralCalculator::calculate(engine, this->scalarBasis().shellSet(), this->scalarBasis().shellSet());


        // For each of the components of the operator, place the scalar basis representations into the spinor basis representation.
        std::array<SquareMatrix<ResultScalar>, N> fs;
        for (size_t i = 0; i < N; i++) {
            fs[i] = SquareMatrix<ResultScalar>(F[i]);
        }


        // 2. Transform using the current coefficient matrix.
        StorageArray<SquareMatrix<ResultScalar>, Vectorizer> array {fs, vectorizer};
        ResultOperator op {array};  // 'op' for 'operator'.
        op.transform(this->expansion());
        return op;
    }


    /**
     *  Quantize the orbital Zeeman operator in this spin-orbital basis.
     * 
     *  @param op               The (first-quantized) orbital Zeeman operator.
     * 
     *  @return The orbital Zeeman operator expressed in this spin-orbital basis.
     */
    template <typename Z = Shell>
    auto quantize(const OrbitalZeemanOperator& op) const -> enable_if_t<std::is_same<Z, LondonGTOShell>::value, USQOneElectronOperatorComponent<product_t<OrbitalZeemanOperator::Scalar, ExpansionScalar>, OrbitalZeemanOperator::Vectorizer>> {

        // Return the orbital Zeeman operator as a contraction beween the magnetic field and the angular momentum operator.
        const auto L = this->quantize(op.angularMomentum());
        const auto& B = op.magneticField().strength();
        return 0.5 * L.dot(B);
    }


    /**
     *  Quantize the diamagnetic operator in this spin-orbital basis.
     * 
     *  @param op               The (first-quantized) diamagnetic operator.
     * 
     *  @return The diamagnetic operator expressed in this spin-orbital basis.
     */
    template <typename Z = Shell>
    auto quantize(const DiamagneticOperator& op) const -> enable_if_t<std::is_same<Z, LondonGTOShell>::value, USQOneElectronOperatorComponent<product_t<DiamagneticOperator::Scalar, ExpansionScalar>, DiamagneticOperator::Vectorizer>> {

        using ResultScalar = product_t<DiamagneticOperator::Scalar, ExpansionScalar>;
        using ResultOperator = USQOneElectronOperatorComponent<ResultScalar, DiamagneticOperator::Vectorizer>;


        // Return the diamagnetic operator as a contraction beween the magnetic field and the electronic quadrupole operator.
        const auto Q = this->quantize(op.electronicQuadrupole()).allParameters();
        const auto& B = op.magneticField().strength();

        // Prepare some variables.
        const auto& B_x = B(CartesianDirection::x);
        const auto& B_y = B(CartesianDirection::y);
        const auto& B_z = B(CartesianDirection::z);

        const auto& Q_xx = Q[DyadicCartesianDirection::xx];
        const auto& Q_xy = Q[DyadicCartesianDirection::xy];
        const auto& Q_xz = Q[DyadicCartesianDirection::xz];
        const auto& Q_yy = Q[DyadicCartesianDirection::yy];
        const auto& Q_yz = Q[DyadicCartesianDirection::yz];
        const auto& Q_zz = Q[DyadicCartesianDirection::zz];


        SquareMatrix<ResultScalar> D_par = 0.125 * ((std::pow(B_y, 2) + std::pow(B_z, 2)) * Q_xx +
                                                    (std::pow(B_x, 2) + std::pow(B_z, 2)) * Q_yy +
                                                    (std::pow(B_x, 2) + std::pow(B_y, 2)) * Q_zz -
                                                    2 * B_x * B_y * Q_xy -
                                                    2 * B_x * B_z * Q_xz -
                                                    2 * B_y * B_z * Q_yz);

        return ResultOperator {D_par};
    }


    /**
     *  MARK: Mulliken domain
     */

    /**
     *  Partition this set of spin-orbitals related to one of the components of an unrestricted spin-orbital basis according to the Mulliken partitioning scheme.
     * 
     *  @param selector             A function that returns true for basis functions that should be included the Mulliken domain.
     * 
     *  @return A `UMullikenDomainComponent` for the AOs selected by the supplied selector function.
     */
    UMullikenDomainComponent<ExpansionScalar> mullikenDomain(const std::function<bool(const BasisFunction&)>& selector) const {

        // FIXME: Try to move this API to SimpleSpinOrbitalBasis.

        const auto ao_indices = this->scalarBasis().basisFunctionIndices(selector);
        return UMullikenDomainComponent<ExpansionScalar> {ao_indices, ao_indices.size()};
    }


    /**
     *  Partition this set of spin-orbitals related to one of the components of an unrestricted spin-orbital basis according to the Mulliken partitioning scheme.
     * 
     *  @param selector             A function that returns true for shells that should be included the Mulliken domain.
     * 
     *  @return A `UMullikenDomainComponent` for the AOs selected by the supplied selector function.
     */
    UMullikenDomainComponent<ExpansionScalar> mullikenDomain(const std::function<bool(const Shell&)>& selector) const {

        // FIXME: Try to move this API to SimpleSpinOrbitalBasis.

        const auto ao_indices = this->scalarBasis().basisFunctionIndices(selector);
        return UMullikenDomainComponent<ExpansionScalar> {ao_indices, ao_indices.size()};
    }
};


/*
 *  MARK: SpinorBasisTraits
 */

/**
 *  A type that provides compile-time information on spinor bases that is otherwise not accessible through a public class alias.
 */
template <typename ExpansionScalar, typename Shell>
struct SpinorBasisTraits<USpinOrbitalBasisComponent<ExpansionScalar, Shell>> {

    // The type of transformation that is naturally related to a `USpinOrbitalBasisComponent`.
    using Transformation = UTransformationComponent<ExpansionScalar>;

    // The second-quantized representation of the overlap operator related to a `USpinOrbitalBasisComponent`.
    using SQOverlapOperator = ScalarUSQOneElectronOperatorComponent<ExpansionScalar>;
};


/*
 *  MARK: BasisTransformableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `BasisTransformable`.
 */
template <typename _ExpansionScalar, typename _Shell>
struct BasisTransformableTraits<USpinOrbitalBasisComponent<_ExpansionScalar, _Shell>> {

    // The type of transformation that is naturally related to a `USpinOrbitalBasisComponent`.
    using Transformation = UTransformationComponent<_ExpansionScalar>;
};


/*
 *  MARK: JacobiRotatableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `JacobiRotatable`.
 */
template <typename _ExpansionScalar, typename _Shell>
struct JacobiRotatableTraits<USpinOrbitalBasisComponent<_ExpansionScalar, _Shell>> {

    // The type of Jacobi rotation that is naturally related to a `USpinOrbitalBasisComponent`.
    using JacobiRotationType = JacobiRotation;
};


}  // namespace GQCP
