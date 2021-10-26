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


#include "Basis/BiorthogonalBasis/GLowdinPairingBasis.hpp"
#include "Basis/BiorthogonalBasis/RLowdinPairingBasis.hpp"
#include "Basis/BiorthogonalBasis/ULowdinPairingBasis.hpp"
#include "Basis/NonOrthogonalBasis/GNonOrthogonalStateBasis.hpp"
#include "Basis/NonOrthogonalBasis/RNonOrthogonalStateBasis.hpp"
#include "Basis/NonOrthogonalBasis/UNonOrthogonalStateBasis.hpp"
#include "DensityMatrix/G1DM.hpp"
#include "DensityMatrix/Orbital1DM.hpp"
#include "DensityMatrix/SpinResolved1DM.hpp"
#include "DensityMatrix/SpinResolved1DMComponent.hpp"
#include "Mathematical/Representation/Matrix.hpp"
#include "Utilities/aliases.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/dynamic_bitset.hpp>


namespace GQCP {


/**
 *  A class that represents an expansion inside a non-orthogonal basis.
 *
 *  @tparam _Scalar                       The scalar type of the expansion coefficients: real or complex.
 *  @tparam _NonOrthogonalBasis           The type of non-orthogonal basis.
 */
template <typename _Scalar, typename _NonOrthogonalBasis>
class NOCIExpansion {
public:
    // The scalar type of the expansion coefficients: real or complex.
    using Scalar = _Scalar;

    // The type of the non-orthogonal basis.
    using NonOrthogonalBasis = _NonOrthogonalBasis;


private:
    // The non-orthpgonal basis with respect to which the coefficients are defined.
    NonOrthogonalBasis non_orthogonal_basis;

    // The expansion coefficients.
    VectorX<Scalar> expansion_coefficients;


public:
    /**
     *  MARK: Constructors
     */

    /**
     *  Construct a NOCI expansion inside the given non-orthogonal basis, with corresponding expansion coefficients.
     *
     *  @param non_orthogonal_basis            The non-orthogonal basis with respect to which the coefficients are defined.
     *  @param coefficients                    The expansion coefficients.
     */
    NOCIExpansion(const NonOrthogonalBasis& non_orthogonal_basis, const VectorX<Scalar>& coefficients) :
        non_orthogonal_basis {non_orthogonal_basis},
        expansion_coefficients {coefficients} {}


    /**
     * MARK: Access
     */

    /**
     *  Access a coefficient of the expansion.
     *
     *  @param i    The index (address) of the coefficient that should be obtained.
     *
     *  @return The i-th expansion coefficient.
     */
    Scalar coefficient(const size_t i) const { return this->expansion_coefficients(i); }

    /**
     *  @return The expansion coefficients of this function model.
     */
    const VectorX<Scalar>& coefficients() const { return this->expansion_coefficients; }

    /**
     *  @return The non-orthogonal basis that is related to this wave function model.
     */
    const NonOrthogonalBasis& nonOrthogonalStateBasis() const { return this->non_orthogonal_basis; }


    /**
     * MARK: Density maytrices for `GNonOrthogonalStateBases`
     */

    /**
     *  Calculate the general one-electron density matrix for a real valued expansion in a generalized non-orthogonal basis.
     *
     *  @return The generalized one-electron density matrix.
     */
    template <typename Z1 = Scalar, typename Z2 = NonOrthogonalBasis>
    enable_if_t<std::is_same<Z1, double>::value && std::is_same<Z2, GNonOrthogonalStateBasis<double>>::value, G1DM<double>> calculate1DM() const {

        // Initialize a zero matrix.
        SquareMatrix<double> D = SquareMatrix<double>::Zero(this->nonOrthogonalStateBasis().basisStateDimension());

        // Loop over all bra-ket combinations in the non-orthogonal basis.
        for (size_t i = 0; i < this->nonOrthogonalStateBasis().numberOfBasisStates(); i++) {

            // Initialize the parameters that are identical for the bra and the ket.
            const auto occupied_orbitals = this->nonOrthogonalStateBasis().numberOfOccupiedOrbitals();
            const auto S = this->nonOrthogonalStateBasis().metric();

            for (size_t j = 0; j < this->nonOrthogonalStateBasis().numberOfBasisStates(); j++) {

                // The first step is to create a biorthogonal basis from the two states that are being looped over.
                const GLowdinPairingBasis<double> lowdin_pairing_basis {this->nonOrthogonalStateBasis().basisState(i), this->nonOrthogonalStateBasis().basisState(j), S, occupied_orbitals};

                // Calculate the transition 1-DM.
                const auto D_xw = lowdin_pairing_basis.transition1DM().matrix();
                const auto total_coefficient = this->coefficient(i) * this->coefficient(j);

                // Calculate the 1DM contribution of the bra and ket, add them to the total matrix and repeat the process.
                D += total_coefficient * D_xw;
            }
        }
        return G1DM<double> {D};
    }


    /**
     *  Calculate the general one-electron density matrix for a complex valued expansion in a generalized non-orthogonal basis.
     *
     *  @return The generalized one-electron density matrix.
     */
    template <typename Z1 = Scalar, typename Z2 = NonOrthogonalBasis>
    enable_if_t<std::is_same<Z1, complex>::value && std::is_same<Z2, GNonOrthogonalStateBasis<complex>>::value, G1DM<complex>> calculate1DM() const {

        // Initialize a zero matrix.
        SquareMatrix<complex> D = SquareMatrix<complex>::Zero(this->nonOrthogonalStateBasis().basisStateDimension());

        // Loop over all bra-ket combinations in the non-orthogonal basis.
        for (size_t i = 0; i < this->nonOrthogonalStateBasis().numberOfBasisStates(); i++) {

            // Initialize the parameters that are identical for the bra and the ket.
            const auto occupied_orbitals = this->nonOrthogonalStateBasis().numberOfOccupiedOrbitals();
            const auto S = this->nonOrthogonalStateBasis().metric();

            for (size_t j = 0; j < this->nonOrthogonalStateBasis().numberOfBasisStates(); j++) {

                // The first step is to create a biorthogonal basis from the two states that are being looped over.
                const GLowdinPairingBasis<complex> lowdin_pairing_basis {this->nonOrthogonalStateBasis().basisState(i), this->nonOrthogonalStateBasis().basisState(j), S, occupied_orbitals};

                // Calculate the transition 1-DM.
                const auto D_xw = lowdin_pairing_basis.transition1DM().matrix();
                const auto total_coefficient = std::conj(this->coefficient(i)) * this->coefficient(j);

                // Calculate the 1DM contribution of the bra and ket, add them to the total matrix and repeat the process.
                D += total_coefficient * D_xw;
            }
        }
        return G1DM<complex> {D};
    }


    /**
     * MARK: Density maytrices for `RNonOrthogonalStateBases`
     */

    /**
     *  Calculate the orbital one-electron density matrix for a real valued expansion in a restricted non-orthogonal basis.
     *
     *  @return The orbital one-electron density matrix.
     */
    template <typename Z1 = Scalar, typename Z2 = NonOrthogonalBasis>
    enable_if_t<std::is_same<Z1, double>::value && std::is_same<Z2, RNonOrthogonalStateBasis<double>>::value, Orbital1DM<double>> calculate1DM() const {

        // Initialize a zero matrix.
        SquareMatrix<double> D = SquareMatrix<double>::Zero(this->nonOrthogonalStateBasis().basisStateDimension());

        // Loop over all bra-ket combinations in the non-orthogonal basis.
        for (size_t i = 0; i < this->nonOrthogonalStateBasis().numberOfBasisStates(); i++) {

            // Initialize the parameters that are identical for the bra and the ket.
            const auto occupied_orbitals = this->nonOrthogonalStateBasis().numberOfOccupiedOrbitals();
            const auto S = this->nonOrthogonalStateBasis().metric();

            for (size_t j = 0; j < this->nonOrthogonalStateBasis().numberOfBasisStates(); j++) {

                // The first step is to create a biorthogonal basis from the two states that are being looped over.
                const RLowdinPairingBasis<double> lowdin_pairing_basis {this->nonOrthogonalStateBasis().basisState(i), this->nonOrthogonalStateBasis().basisState(j), S, occupied_orbitals};

                // Calculate the transition 1-DM.
                const auto D_xw = lowdin_pairing_basis.transition1DM().matrix();
                const auto total_coefficient = this->coefficient(i) * this->coefficient(j);

                // Calculate the 1DM contribution of the bra and ket, add them to the total matrix and repeat the process.
                D += total_coefficient * D_xw;
            }
        }
        return Orbital1DM<double> {D};
    }


    /**
     *  Calculate the orbital one-electron density matrix for a complex valued expansion in a restricted non-orthogonal basis.
     *
     *  @return The orbital one-electron density matrix.
     */
    template <typename Z1 = Scalar, typename Z2 = NonOrthogonalBasis>
    enable_if_t<std::is_same<Z1, complex>::value && std::is_same<Z2, RNonOrthogonalStateBasis<complex>>::value, Orbital1DM<complex>> calculate1DM() const {

        // Initialize a zero matrix.
        SquareMatrix<complex> D = SquareMatrix<complex>::Zero(this->nonOrthogonalStateBasis().basisStateDimension());

        // Loop over all bra-ket combinations in the non-orthogonal basis.
        for (size_t i = 0; i < this->nonOrthogonalStateBasis().numberOfBasisStates(); i++) {

            // Initialize the parameters that are identical for the bra and the ket.
            const auto occupied_orbitals = this->nonOrthogonalStateBasis().numberOfOccupiedOrbitals();
            const auto S = this->nonOrthogonalStateBasis().metric();

            for (size_t j = 0; j < this->nonOrthogonalStateBasis().numberOfBasisStates(); j++) {

                // The first step is to create a biorthogonal basis from the two states that are being looped over.
                const RLowdinPairingBasis<complex> lowdin_pairing_basis {this->nonOrthogonalStateBasis().basisState(i), this->nonOrthogonalStateBasis().basisState(j), S, occupied_orbitals};

                // Calculate the transition 1-DM.
                const auto D_xw = lowdin_pairing_basis.transition1DM().matrix();
                const auto total_coefficient = std::conj(this->coefficient(i)) * this->coefficient(j);

                // Calculate the 1DM contribution of the bra and ket, add them to the total matrix and repeat the process.
                D += total_coefficient * D_xw;
            }
        }
        return Orbital1DM<complex> {D};
    }


    /**
     * MARK: Density maytrices for `UNonOrthogonalStateBases`
     */

    /**
     *  Calculate the spin-resolved one-electron density matrix for a real valued expansion in a unrestricted non-orthogonal basis.
     *
     *  @return The spin-resolved one-electron density matrix.
     */
    template <typename Z1 = Scalar, typename Z2 = NonOrthogonalBasis>
    enable_if_t<std::is_same<Z1, double>::value && std::is_same<Z2, UNonOrthogonalStateBasis<double>>::value, SpinResolved1DM<double>> calculate1DM() const {

        // Initialize a zero matrix.
        SquareMatrix<double> D_a = SquareMatrix<double>::Zero(this->nonOrthogonalStateBasis().basisStateDimension());
        SquareMatrix<double> D_b = SquareMatrix<double>::Zero(this->nonOrthogonalStateBasis().basisStateDimension());

        // Loop over all bra-ket combinations in the non-orthogonal basis.
        for (size_t i = 0; i < this->nonOrthogonalStateBasis().numberOfBasisStates(); i++) {

            // Initialize the parameters that are identical for the bra and the ket.
            const auto occupied_alpha_orbitals = this->nonOrthogonalStateBasis().numberOfOccupiedOrbitals().alpha();
            const auto occupied_beta_orbitals = this->nonOrthogonalStateBasis().numberOfOccupiedOrbitals().beta();
            const auto S = this->nonOrthogonalStateBasis().metric();

            for (size_t j = 0; j < this->nonOrthogonalStateBasis().numberOfBasisStates(); j++) {

                // The first step is to create a biorthogonal basis from the two states that are being looped over.
                const ULowdinPairingBasis<double> lowdin_pairing_basis {this->nonOrthogonalStateBasis().basisState(i), this->nonOrthogonalStateBasis().basisState(j), S, occupied_alpha_orbitals, occupied_beta_orbitals};

                // Calculate the transition 1-DM.
                const auto D_xw_a = lowdin_pairing_basis.transition1DM().alpha().matrix();
                const auto D_xw_b = lowdin_pairing_basis.transition1DM().beta().matrix();
                const auto total_coefficient = this->coefficient(i) * this->coefficient(j);

                // Calculate the 1DM contribution of the bra and ket, add them to the total matrix and repeat the process.
                D_a += total_coefficient * D_xw_a;
                D_b += total_coefficient * D_xw_b;
            }
        }
        return SpinResolved1DM<double> {D_a, D_b};
    }


    /**
     *  Calculate the spin-resolved one-electron density matrix for a complex valued expansion in a unrestricted non-orthogonal basis.
     *
     *  @return The spin-resolved one-electron density matrix.
     */
    template <typename Z1 = Scalar, typename Z2 = NonOrthogonalBasis>
    enable_if_t<std::is_same<Z1, complex>::value && std::is_same<Z2, UNonOrthogonalStateBasis<complex>>::value, SpinResolved1DM<complex>> calculate1DM() const {

        // Initialize a zero matrix.
        SquareMatrix<complex> D_a = SquareMatrix<complex>::Zero(this->nonOrthogonalStateBasis().basisStateDimension());
        SquareMatrix<complex> D_b = SquareMatrix<complex>::Zero(this->nonOrthogonalStateBasis().basisStateDimension());

        // Loop over all bra-ket combinations in the non-orthogonal basis.
        for (size_t i = 0; i < this->nonOrthogonalStateBasis().numberOfBasisStates(); i++) {

            // Initialize the parameters that are identical for the bra and the ket.
            const auto occupied_alpha_orbitals = this->nonOrthogonalStateBasis().numberOfOccupiedOrbitals().alpha();
            const auto occupied_beta_orbitals = this->nonOrthogonalStateBasis().numberOfOccupiedOrbitals().beta();
            const auto S = this->nonOrthogonalStateBasis().metric();

            for (size_t j = 0; j < this->nonOrthogonalStateBasis().numberOfBasisStates(); j++) {

                // The first step is to create a biorthogonal basis from the two states that are being looped over.
                const ULowdinPairingBasis<complex> lowdin_pairing_basis {this->nonOrthogonalStateBasis().basisState(i), this->nonOrthogonalStateBasis().basisState(j), S, occupied_alpha_orbitals, occupied_beta_orbitals};

                // Calculate the transition 1-DM.
                const auto D_xw_a = lowdin_pairing_basis.transition1DM().alpha().matrix();
                const auto D_xw_b = lowdin_pairing_basis.transition1DM().beta().matrix();
                const auto total_coefficient = this->coefficient(i) * this->coefficient(j);


                // Calculate the 1DM contribution of the bra and ket, add them to the total matrix and repeat the process.
                D_a += total_coefficient * D_xw_a;
                D_b += total_coefficient * D_xw_b;
            }
        }
        return SpinResolved1DM<complex> {D_a, D_b};
    }
};


}  // namespace GQCP
