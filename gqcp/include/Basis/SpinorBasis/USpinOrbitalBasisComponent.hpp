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


#include "Basis/MullikenPartitioning/UMullikenPartitioningComponent.hpp"
#include "Basis/SpinorBasis/SimpleSpinOrbitalBasis.hpp"
#include "Basis/Transformations/JacobiRotation.hpp"
#include "Basis/Transformations/UTransformationComponent.hpp"
#include "Operator/SecondQuantized/USQOneElectronOperatorComponent.hpp"


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
     *  MARK: Quantization of first-quantized operators
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
    template <typename FQOneElectronOperator>
    auto quantize(const FQOneElectronOperator& fq_one_op) const -> USQOneElectronOperatorComponent<product_t<typename FQOneElectronOperator::Scalar, ExpansionScalar>, typename FQOneElectronOperator::Vectorizer> {

        // FIXME: Try to move this API to SimpleSpinOrbitalBasis.
        using ResultScalar = product_t<typename FQOneElectronOperator::Scalar, ExpansionScalar>;
        using ResultOperator = USQOneElectronOperatorComponent<ResultScalar, typename FQOneElectronOperator::Vectorizer>;

        const auto f_par = IntegralCalculator::calculateLibintIntegrals(fq_one_op, this->scalarBasis());  // Currently f is expressed in the AO basis. 'par' for 'parameters'

        ResultOperator f {f_par};
        f.transform(this->expansion());  // Now, f is expressed in the spin-orbital basis.
        return f;
    }


    /**
     *  MARK: Mulliken partitioning
     */

    /**
     *  Partition this set of spin-orbitals related to one of the components of an unrestricted spin-orbital basis according to the Mulliken partitioning scheme.
     * 
     *  @param selector             A function that returns true for basis functions that should be included the Mulliken partitioning.
     * 
     *  @return A `UMullikenPartitioningComponent` for the AOs selected by the supplied selector function.
     */
    UMullikenPartitioningComponent<ExpansionScalar> mullikenPartitioning(const std::function<bool(const BasisFunction&)>& selector) const {

        // FIXME: Try to move this API to SimpleSpinOrbitalBasis.

        const auto ao_indices = this->scalarBasis().basisFunctionIndices(selector);
        return UMullikenPartitioningComponent<ExpansionScalar> {ao_indices, this->expansion()};
    }


    /**
     *  Partition this set of spin-orbitals related to one of the components of an unrestricted spin-orbital basis according to the Mulliken partitioning scheme.
     * 
     *  @param selector             A function that returns true for shells that should be included the Mulliken partitioning.
     * 
     *  @return A `UMullikenPartitioningComponent` for the AOs selected by the supplied selector function.
     */
    UMullikenPartitioningComponent<ExpansionScalar> mullikenPartitioning(const std::function<bool(const Shell&)>& selector) const {

        // FIXME: Try to move this API to SimpleSpinOrbitalBasis.

        const auto ao_indices = this->scalarBasis().basisFunctionIndices(selector);
        return UMullikenPartitioningComponent<ExpansionScalar> {ao_indices, this->expansion()};
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
