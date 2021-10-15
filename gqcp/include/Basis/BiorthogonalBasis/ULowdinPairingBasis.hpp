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


#include "Basis/BiorthogonalBasis/ULowdinPairingBasis.hpp"
#include "Basis/Transformations/UTransformation.hpp"
#include "Operator/SecondQuantized/USQOneElectronOperator.hpp"
#include "QuantumChemical/SpinResolved.hpp"
#include "QuantumChemical/SpinResolvedBase.hpp"


namespace GQCP {


/*
 *  MARK: ULowdinPairingBasis implementation
 */

/**
 *  A type that encapsulates biorthogonal Lowdin pairing basis, split in its alpha and beta components.
 *
 *  @tparam _Scalar         The scalar type used for the expansion coefficients: real or complex.
 */
template <typename _Scalar>
class ULowdinPairingBasis:
    public SpinResolvedBase<ULowdinPairingBasisComponent<_Scalar>, ULowdinPairingBasis<_Scalar>> {
public:
    // The scalar type used for the expansion coefficients: real or complex.
    using Scalar = _Scalar;

    // The type of 'this'.
    using Self = ULowdinPairingBasis<Scalar>;

    // The type component this spin resolved object is made of.
    using ComponentType = typename SpinResolvedBase<ULowdinPairingBasisComponent<Scalar>, Self>::Of;


public:
    /*
     *  MARK: Constructors
     */

    // Inherit `SpinResolvedBase`'s constructors.
    using SpinResolvedBase<ULowdinPairingBasisComponent<Scalar>, ULowdinPairingBasis<Scalar>>::SpinResolvedBase;


    /*
     *  MARK: Constructors
     */

    /**
     *  Create a `ULowdinPairingBasis` from two non-orthogonal states.
     *
     *  @param C_bra                                  The transformation that represents the expansion coefficients of the bra non-orthogonal state. This is a UTransformation in this case.
     *  @param C_ket                                  The transformation that represents the expansion coefficients of the ket non-orthogonal state. This is a UTransformation in this case.
     *  @param S_AO                                   The overlap operator in AO basis, constructed from the `USpinOrbitalBasis` used to calculate the non-orthogonal states.
     *  @param number_of_occupied_alpha_orbitals      The total number of occupied orbitals in the system.
     *  @param number_of_occupied_beta_orbitals       The total number of occupied orbitals in the system.
     *  @param threshold                              The threshold at which a value is verified to be zero or not. Default is 1e-8.
     */
    ULowdinPairingBasis<Scalar>(const UTransformation<Scalar>& C_bra, const UTransformation<Scalar>& C_ket, const USQOneElectronOperator<Scalar>& S_AO, const size_t number_of_occupied_alpha_orbitals, const size_t number_of_occupied_beta_orbitals, const double threshold = 1e-8) :
        ULowdinPairingBasis(ULowdinPairingBasisComponent<Scalar> {C_bra.alpha(), C_ket.alpha(), S_AO.alpha(), number_of_occupied_alpha_orbitals, threshold},
                            ULowdinPairingBasisComponent<Scalar> {C_bra.beta(), C_ket.beta(), S_AO.beta(), number_of_occupied_beta_orbitals, threshold}) {};
};


}  // namespace GQCP
