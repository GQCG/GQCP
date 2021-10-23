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


#include "Basis/BiorthogonalBasis/ULowdinPairingBasisComponent.hpp"
#include "Basis/Transformations/UTransformation.hpp"
#include "DensityMatrix/SpinResolved1DM.hpp"
#include "Operator/SecondQuantized/USQOneElectronOperator.hpp"
#include "QuantumChemical/Spin.hpp"
#include "QuantumChemical/SpinResolved.hpp"
#include "QuantumChemical/SpinResolvedBase.hpp"


namespace GQCP {


/*
 *  MARK: ULowdinPairingBasis implementation
 */

/**
 *  A type that encapsulates biorthogonal Löwdin pairing basis, split in its alpha and beta components.
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

    // The type of matrix naturally associated with the components of a `ULowdinPairingBasis`.
    using Matrix = MatrixX<Scalar>;

    // The vectors associated with the scalar of the expansion coefficients.
    using Vector = VectorX<Scalar>;


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
     *  @param number_of_occupied_alpha_orbitals      The total number of occupied alpha orbitals in the system.
     *  @param number_of_occupied_beta_orbitals       The total number of occupied beta orbitals in the system.
     *  @param threshold                              The threshold at which a value is verified to be zero or not. The default is 1e-8.
     */
    ULowdinPairingBasis<Scalar>(const UTransformation<Scalar>& C_bra, const UTransformation<Scalar>& C_ket, const ScalarUSQOneElectronOperator<Scalar>& S_AO, const size_t number_of_occupied_alpha_orbitals, const size_t number_of_occupied_beta_orbitals, const double threshold = 1e-8) :
        ULowdinPairingBasis(ULowdinPairingBasisComponent<Scalar> {C_bra.alpha(), C_ket.alpha(), S_AO.alpha(), number_of_occupied_alpha_orbitals, threshold},
                            ULowdinPairingBasisComponent<Scalar> {C_bra.beta(), C_ket.beta(), S_AO.beta(), number_of_occupied_beta_orbitals, threshold}) {};


    /*
     *  MARK: Access
     */

    /**
     * Return the biorthogonalized expansion coefficients of both the bra and the ket as a pair.
     *
     * @return The biorthogonalized expansion coefficients.
     *
     * @note   Only the occupied expansion coefficients are biorthogonalized, so only these coefficients are returned.
     */
    const SpinResolved<std::pair<Matrix, Matrix>> biorthogonalExpansion() const { return SpinResolved<std::pair<Matrix, Matrix>> {this->alpha().biorthogonalExpansion(), this->beta().biorthogonalExpansion()}; }


    /**
     * Return the biorthogonalized expansion coefficients of both the bra and the ket as a pair.
     *
     * @param sigma     The spin sigma component: alpha or beta.
     *
     * @return The biorthogonalized expansion coefficients of the alpha or beta component.
     *
     * @note   Only the occupied expansion coefficients are biorthogonalized, so only these coefficients are returned.
     */
    const std::pair<Matrix, Matrix>& biorthogonalExpansion(const Spin sigma) const { return this->component(sigma).biorthogonalExpansion(); }


    /**
     * Return the biorthogonalized expansion coefficients of the bra.
     *
     * @return The biorthogonalized expansion coefficients belonging to the bra.
     *
     * @note   Only the occupied expansion coefficients are biorthogonalized, so only these coefficients are returned.
     */
    const SpinResolved<Matrix>& biorthogonalBraExpansion() const { return SpinResolved<Matrix> {this->alpha().biorthogonalExpansion().first, this->beta().biorthogonalExpansion().first}; }


    /**
     * Return the biorthogonalized expansion coefficients of the bra.
     *
     * @param sigma     The spin sigma component: alpha or beta.
     *
     * @return The biorthogonalized expansion coefficients belonging to the bra of the alpha or beta component.
     *
     * @note   Only the occupied expansion coefficients are biorthogonalized, so only these coefficients are returned.
     */
    const Matrix& biorthogonalBraExpansion(const Spin sigma) const { return this->biorthogonalExpansion(sigma).first; }


    /**
     * Return the biorthogonalized expansion coefficients of the bra.
     *
     * @return The biorthogonalized expansion coefficients belonging to the bra.
     *
     * @note   Only the occupied expansion coefficients are biorthogonalized, so only these coefficients are returned.
     */
    const SpinResolved<Matrix>& biorthogonalKetExpansion() const { return SpinResolved<Matrix> {this->alpha().biorthogonalExpansion().second, this->beta().biorthogonalExpansion().second}; }


    /**
     * Return the biorthogonalized expansion coefficients of the ket.
     *
     * @param sigma     The spin sigma component: alpha or beta.
     *
     * @return The biorthogonalized expansion coefficients belonging to the ket of the alpha or beta component.
     *
     * @note   Only the occupied expansion coefficients are biorthogonalized, so only these coefficients are returned.
     */
    const Matrix& biorthogonalKetExpansion(const Spin sigma) const { return this->biorthogonalExpansion(sigma).second; }


    /**
     * Return the overlap values of the biorthogonalized expansion coefficients.
     *
     * @return The overlaps of the biorthogonal coefficients.
     */
    const SpinResolved<Vector>& biorthogonalOverlaps() const { return SpinResolved<Vector> {this->alpha().biorthogonalOverlaps(), this->beta().biorthogonalOverlaps()}; }


    /**
     * Return the overlap values of the biorthogonalized expansion coefficients.
     *
     * @param sigma     The spin sigma component: alpha or beta.
     *
     * @return The overlaps of the biorthogonal coefficients of the alpha or beta component.
     */
    const Vector& biorthogonalOverlaps(const Spin sigma) const { return this->component(sigma).biorthogonalOverlaps(); }


    /*
     *  MARK: Overlap
     */

    /**
     * Determine the number of zero overlaps in the biorthogonal overlap vector of a certain spin component.
     *
     * @param sigma     The spin sigma component: alpha or beta.
     *
     * @return The number of zero overlaps in the spin sigma overlap values.
     */
    int numberOfZeroOverlaps(const Spin sigma) const { return this->component(sigma).numberOfZeroOverlaps(); }


    /**
     * Determine the total number of zero overlaps in the biorthogonal overlap vectors.
     *
     * @return The number of zero overlaps in both the spin sigma overlap vectors.
     */
    int numberOfZeroOverlaps() const { return this->numberOfZeroOverlaps(Spin::alpha) + this->numberOfZeroOverlaps(Spin::beta); }

    /**
     * Calculate the reduced overlap. I.e. The biorthogonal overlaps with the zero values get removed and the total remaining overlap is calculated.
     *
     * @param sigma     The spin sigma component: alpha or beta.
     *
     * @return The reduced overlap of the spin sigma component.
     */
    Scalar reducedOverlap(const Spin sigma) const { return this->component(sigma).reducedOverlap(); }


    /**
     * Calculate the total reduced overlap. I.e. The biorthogonal overlaps with the zero values get removed and the total remaining overlap is calculated.
     *
     * @return The total reduced overlap.
     */
    Scalar reducedOverlap() const { return this->reducedOverlap(Spin::alpha) * this->reducedOverlap(Spin::beta); }


    /**
     * Calculate and return the total overlap value of the biorthogonal coefficients of the spin sigma component.
     *
     * @param sigma     The spin sigma component: alpha or beta.
     *
     * @return The total overlap value of the spin sigma component.
     */
    Scalar totalOverlap(const Spin sigma) const { return this->biorthogonalOverlaps(sigma).prod(); }


    /**
     * Calculate and return the total overlap value of the biorthogonal coefficients.
     *
     * @return The total overlap value.
     */
    Scalar totalOverlap() const { return this->totalOverlap(Spin::alpha) * this->totalOverlap(Spin::beta); }


    /**
     * Determine the indices of the zero overlap values in the biorthogonal overlap vector of the spin sigma component.
     *
     * @param sigma     The spin sigma component: alpha or beta.
     *
     * @return The indices of the zero overlap values of the spin sigma component.
     */
    std::vector<size_t> zeroOverlapIndices(const Spin sigma) const { return this->component(sigma).zeroOverlapIndices(); }


    /**
     * Determine the indices of the zero overlap values in the biorthogonal overlap vector of each spin component.
     *
     * @return The indices of the zero overlap values of the spin sigma component.
     */
    std::vector<std::pair<size_t, Spin>> zeroOverlapIndices() const {

        // We merge the alpha and beta zero-index vectors, but we keep track of the spins of each one.
        std::vector<std::pair<size_t, Spin>> pair_vector {};

        for (const auto& index : this->zeroOverlapIndices(Spin::alpha)) {
            pair_vector.push_back(std::pair<size_t, Spin> {index, Spin::alpha});
        }

        for (const auto& index : this->zeroOverlapIndices(Spin::beta)) {
            pair_vector.push_back(std::pair<size_t, Spin> {index, Spin::beta});
        }

        return pair_vector;
    }


    /*
     *  MARK: Density matrices
     */

    /**
     * Calculate and return the co-density matrix for the given occupied orbital index.
     *
     * @param k     The occupied orbital index.
     *
     * @return The co-density matrix at occupied orbital index k.
     *
     * @note This implementation is based on equation 38b from the 2021 paper by Hugh Burton (https://aip.scitation.org/doi/abs/10.1063/5.0045442).
     */
    SpinResolved1DM<Scalar> coDensity(const int k) const { return SpinResolved1DM<Scalar> {SpinResolved1DMComponent<Scalar> {this->alpha().coDensity(k).matrix()}, SpinResolved1DMComponent<Scalar> {this->beta().coDensity(k).matrix()}}; }


    /**
     * Calculate and return the sum of the zero overlap co-density matrix, the weighted co-density matrix and the zero overlap co-density matrix of the biorthogonalized basis of the bra with itself.
     *
     * @return The sum of the co-density matrices.
     *
     * @note This implementation is based on equation 38d from the 2021 paper by Hugh Burton (https://aip.scitation.org/doi/abs/10.1063/5.0045442).
     */
    SpinResolved1DM<Scalar> coDensitySum() const { return SpinResolved1DM<Scalar> {SpinResolved1DMComponent<Scalar> {this->alpha().coDensitySum().matrix()}, SpinResolved1DMComponent<Scalar> {this->beta().coDensitySum().matrix()}}; }

    /**
     * Calculate and return the transition one-electron density matrix. The contractions vary depending on the number of zero overlaps (Burton equation 46 and 47).
     * These matrices are then summed and returned.
     *
     * @return The transition one-electron density matrix.
     *
     * @note This implementation is based on equation 45 from the 2021 paper by Hugh Burton (https://aip.scitation.org/doi/abs/10.1063/5.0045442).
     */
    SpinResolved1DM<Scalar> transition1DM() const { return SpinResolved1DM<Scalar> {SpinResolved1DMComponent<Scalar> {this->alpha().transition1DM().matrix()}, SpinResolved1DMComponent<Scalar> {this->beta().transition1DM().matrix()}}; }


    /**
     * Calculate and return the weighted co-density matrix. It calculates the co-density matrix contributions corresponding to the non-zero overlap values and divides them by that overlap value.
     * These matrices are then summed and returned.
     *
     * @return The weighted co-density matrix.
     *
     * @note This implementation is based on equation 38c from the 2021 paper by Hugh Burton (https://aip.scitation.org/doi/abs/10.1063/5.0045442).
     */
    SpinResolved1DM<Scalar> weightedCoDensity() const { return SpinResolved1DM<Scalar> {SpinResolved1DMComponent<Scalar> {this->alpha().weightedCoDensity().matrix()}, SpinResolved1DMComponent<Scalar> {this->beta().weightedCoDensity().matrix()}}; }


    /**
     * Calculate and return the zero-overlap co-density matrix. It takes the indices of the zero overlaps and calculates the co-density matrix at each of these indices. The sum of those co-density matrices is returned.
     *
     * @return The zero overlap co-density matrix.
     *
     * @note This implementation is based on equation 38a from the 2021 paper by Hugh Burton (https://aip.scitation.org/doi/abs/10.1063/5.0045442).
     */
    SpinResolved1DM<Scalar> zeroOverlapCoDensity() const { return SpinResolved1DM<Scalar> {SpinResolved1DMComponent<Scalar> {this->alpha().zeroOverlapCoDensity().matrix()}, SpinResolved1DMComponent<Scalar> {this->beta().zeroOverlapCoDensity().matrix()}}; }
};


}  // namespace GQCP
