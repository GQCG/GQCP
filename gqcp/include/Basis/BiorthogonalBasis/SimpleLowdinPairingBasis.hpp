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


#include "Mathematical/Representation/Tensor.hpp"
#include "Utilities/CRTP.hpp"
#include "Utilities/Eigen.hpp"
#include "Utilities/complex.hpp"

#include <Eigen/SVD>

#include <numeric>


namespace GQCP {

/*
 *  MARK: LowdinPairingBasisTraits
 */

/**
 *  A type that provides compile-time information on Löwdin pairing bases that is otherwise not accessible through a public class alias.
 */
template <typename LowdinPairingBasis>
struct LowdinPairingBasisTraits {};

/*
 *  MARK: SimpleLowdinPairingBasis
 */

/**
 *  A Löwdin pairing basis formed from two `R/U/GTransformation`s. Two given sets of expansion coefficients are bi-orthogonalised.
 *
 *  @tparam _ExpansionScalar        The scalar type used to represent the expansion coefficients of the given non-orthogonal states: real or complex.
 */
template <typename _Scalar, typename _DerivedLowdinPairingBasis>
class SimpleLowdinPairingBasis {
public:
    // The scalar type used to represent the expansion coefficients of the given non-orthogonal states: real or complex.
    using Scalar = _Scalar;

    // The vectors associated with the scalar of the expansion coefficients.
    using Vector = VectorX<Scalar>;

    // The matrices associated with the scalar of the expansion coefficients.
    using Matrix = MatrixX<Scalar>;

    // The type of the derived Lowdin pairing Basis from this parent class.
    using DerivedLowdinPairingBasis = _DerivedLowdinPairingBasis;

    // The type of transformation that is naturally related to a `LowdinPairingBasis`. Can be R/U/GTransformation.
    using Transformation = typename LowdinPairingBasisTraits<DerivedLowdinPairingBasis>::Transformation;

    // The second-quantized representation of the overlap operator related to the final Löwdin pairing basis.
    using SQOverlapOperator = typename LowdinPairingBasisTraits<DerivedLowdinPairingBasis>::SQOverlapOperator;

    // The density matrix associated with a lowdin paring basis.
    using DM = typename LowdinPairingBasisTraits<DerivedLowdinPairingBasis>::DM;

    // The representation of this `SimpleLowdinPairingBasis`.
    using Self = SimpleLowdinPairingBasis<Scalar, DerivedLowdinPairingBasis>;


protected:
    // The total number of occupied orbitals.
    size_t N;

    // The total number of orbitals.
    size_t M;

    // The threshold used to determine zero values.
    double zero_threshold;

    // The expansion coefficients of the two non orthogonal states. The element `first`corresponds to the bra expansion, the element `second` corresponds to the ket expansion.
    std::pair<Transformation, Transformation> non_orthogonal_state_expansions;

    // The expansion coefficients of the two biorthogonalized states. The element `first`corresponds to the biorthogonal bra expansion, the element `second` corresponds to the biorthogonal ket expansion.
    // The expansion coefficients of the occupied biorthogonal states are saved as raw matrices instead of the corresponding `Transformation`s because they are not square.
    std::pair<Matrix, Matrix> occupied_biorthogonal_state_expansions;

    // The overlap operator in AO basis, constructed from the spinor/spin-orbital basis used to calculate the non-orthogonal states.
    SQOverlapOperator overlap_operator_AO;

    // The vector containing the biorthogonal overlaps.
    Vector biorthogonal_overlaps;

public:
    /*
     *  MARK: Constructors
     */

    /**
     *  Create a `SimpleLowdinPairingBasis` from two non-orthogonal states.
     *
     *  @param C_bra                            The transformation that represents the expansion coefficients of the bra non-orthogonal state.
     *  @param C_ket                            The transformation that represents the expansion coefficients of the ket non-orthogonal state.
     *  @param S_AO                             The overlap operator in AO basis, constructed from the spinor/spin-orbital used to calculate the non-orthogonal states.
     *  @param number_of_occupied_orbitals      The total number of occupied orbitals in the system.
     *  @param threshold                        The threshold at which a value is verified to be zero or not. The default is 1e-8.
     */
    SimpleLowdinPairingBasis<Scalar, DerivedLowdinPairingBasis>(const Transformation& C_bra, const Transformation& C_ket, const SQOverlapOperator& S_AO, const size_t number_of_occupied_orbitals, const double threshold = 1e-8) :
        non_orthogonal_state_expansions {std::pair<Transformation, Transformation> {C_bra, C_ket}},
        overlap_operator_AO {S_AO},
        zero_threshold {threshold},
        N {number_of_occupied_orbitals},
        M {0},
        biorthogonal_overlaps {Vector {0}},
        occupied_biorthogonal_state_expansions {std::pair<Matrix, Matrix> {Matrix::Zero(1, 1), Matrix::Zero(1, 1)}} {

        // Initialize the number of orbitals.
        this->M = this->overlap_operator_AO.parameters().rows();

        // To determine the biorthogonal overlaps, we only need the coefficients of the occupied orbitals.
        Matrix C_bra_occupied {this->M, this->N};
        Matrix C_ket_occupied {this->M, this->N};

        for (size_t i = 0; i < this->M; i++) {
            for (size_t j = 0; j < this->N; j++) {
                C_bra_occupied(i, j) = C_bra.matrix()(i, j);
                C_ket_occupied(i, j) = C_ket.matrix()(i, j);
            }
        }
        // std::cout << "-----bra-----" << std::endl;
        // std::cout << C_bra_occupied << std::endl;
        // std::cout << "-----ket-----" << std::endl;
        // std::cout << C_ket_occupied << std::endl;

        // Now that we have the occupied expansions of the bra and the ket, we calculate their overlap.
        const auto occupied_orbital_overlap = C_bra_occupied.transpose().conjugate() * S_AO.parameters() * C_ket_occupied;

        // std::cout << "-----S-bra-ket-----" << std::endl;
        // std::cout << occupied_orbital_overlap << std::endl;

        // Perform a singular value decomposition (SVD) on the occupied orbital overlap, in order to gain the biorthogonal overlaps.
        // We require the full matrices from the SVD, not the so-called `thin` matrices.
        Eigen::JacobiSVD<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>> svd(occupied_orbital_overlap, Eigen::ComputeFullU | Eigen::ComputeFullV);

        // We need a two dimensional tensor representation of the provided expansions in order to perform a contraction later on.
        Eigen::TensorMap<Eigen::Tensor<const Scalar, 2>> occ_bra_map {C_bra_occupied.data(), C_bra_occupied.rows(), C_bra_occupied.cols()};
        Tensor<Scalar, 2> C_bra_occupied_tensor = Tensor<Scalar, 2>(occ_bra_map);

        Eigen::TensorMap<Eigen::Tensor<const Scalar, 2>> occ_ket_map {C_ket_occupied.data(), C_ket_occupied.rows(), C_ket_occupied.cols()};
        Tensor<Scalar, 2> C_ket_occupied_tensor = Tensor<Scalar, 2>(occ_ket_map);

        // std::cout << "-----U-----" << std::endl;
        // std::cout << svd.matrixU() << std::endl;
        // std::cout << "-----V-----" << std::endl;
        // std::cout << svd.matrixV() << std::endl;

        // Perform the contractions in order to biorthogonalize the occupied bra and ket expansions.
        // The SVD, in contrast to numpy, returns the matrices U and V. No adjoint or transpose is necessary.
        Matrix biorthogonal_bra_occupied = C_bra_occupied_tensor.template einsum<1>("ui,ij->uj", svd.matrixU()).asMatrix();
        Matrix biorthogonal_ket_occupied = C_ket_occupied_tensor.template einsum<1>("ui,ij->uj", svd.matrixV()).asMatrix();
        // std::cout << "-----biorth-bra-no-phase-----" << std::endl;
        // std::cout << biorthogonal_bra_occupied << std::endl;
        // std::cout << "-----biorth-ket-no-phase-----" << std::endl;
        // std::cout << biorthogonal_ket_occupied << std::endl;
        // Correct the biorthogonal expansions for the phase factor.
        biorthogonal_bra_occupied.col(0) = biorthogonal_bra_occupied.col(0) * svd.matrixU().transpose().determinant();
        biorthogonal_ket_occupied.col(0) = biorthogonal_ket_occupied.col(0) * svd.matrixV().transpose().determinant();
        // std::cout << "-----det(U.T)-----" << std::endl;
        // std::cout << svd.matrixU().transpose().determinant() << std::endl;
        // std::cout << "-----det(V.T)-----" << std::endl;
        // std::cout << svd.matrixV().transpose().determinant() << std::endl;
        // We now have the biorthogonal expansion coefficients.
        this->occupied_biorthogonal_state_expansions = std::pair<Matrix, Matrix> {biorthogonal_bra_occupied, biorthogonal_ket_occupied};

        // We will perform another check to see whether the biorthogonalization procedure was correct.
        // The overlap between the biorthogonal states should be the same when calculated by multiplying the biorthogonal overlaps as when taking the determinant of the overlap matrix X of the occupied orbitals only.
        const auto X = this->occupied_biorthogonal_state_expansions.first.conjugate().transpose() * this->overlap_operator_AO.parameters() * this->occupied_biorthogonal_state_expansions.second;

        // Generate the occupied only overlap matrix X_occupied.
        Matrix X_occupied {this->N, this->N};
        for (size_t i = 0; i < this->N; i++) {
            for (size_t j = 0; j < this->N; j++) {
                X_occupied(i, j) = X(i, j);
            }
        }

        // The singular values are the biorthogonal overlaps.
        const auto new_overlaps = X.diagonal();
        this->biorthogonal_overlaps = new_overlaps;

        // Determine the procduct of the biorthogonal overlaps.
        const auto overlap = this->biorthogonal_overlaps.prod();

        // We now check the aforementioned condition.
        if (std::abs(overlap - X_occupied.determinant()) > 1e-8) {
            throw std::invalid_argument("LowdinPairingBasis<Scalar>(const Transformation& C_bra, const Transformation& C_ket, const ScalarGSQOneElectronOperator<Scalar>& S_AO, const size_t number_of_electrons, const double threshold): The given parameters lead to a wrong overlap calculation.");
        }
    }


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
    const std::pair<Matrix, Matrix>& biorthogonalExpansion() const { return this->occupied_biorthogonal_state_expansions; }

    /**
     * Return the biorthogonalized expansion coefficients of the bra.
     *
     * @return The biorthogonalized expansion coefficients belonging to the bra.
     *
     * @note   Only the occupied expansion coefficients are biorthogonalized, so only these coefficients are returned.
     */
    const Matrix& biorthogonalBraExpansion() const { return this->biorthogonalExpansion().first; }

    /**
     * Return the biorthogonalized expansion coefficients of the ket.
     *
     * @return The biorthogonalized expansion coefficients belonging to the ket.
     *
     * @note   Only the occupied expansion coefficients are biorthogonalized, so only these coefficients are returned.
     */
    const Matrix& biorthogonalKetExpansion() const { return this->biorthogonalExpansion().second; }

    /**
     * Return the overlap values of the biorthogonalized expansion coefficients.
     *
     * @return The overlap values of the biorthogonal coefficients.
     */
    const Vector& biorthogonalOverlaps() const { return this->biorthogonal_overlaps; }

    /**
     * Return the threshold used to compare values to zero associated with this biorthogonal Löwdin Pairing basis.
     *
     * @return The overlap values of the biorthogonal coefficients.
     */
    const double& threshold() const { return this->zero_threshold; }


    /*
     *  MARK: Overlap
     */

    /**
     * Determine the number of zero overlaps in the biorthogonal overlap vector.
     *
     * @return The number of zero overlaps.
     */
    template <typename Z = Scalar>
    enable_if_t<std::is_same<Z, double>::value, int> numberOfZeroOverlaps() const {

        // The count starets at zero.
        int number_of_zeros = 0;

        // Check all overlap values and increase the count if the overlap value is zero.
        for (int i = 0; i < this->biorthogonalOverlaps().rows(); i++) {
            if (std::abs(this->biorthogonalOverlaps()[i]) < this->zero_threshold) {
                number_of_zeros += 1;
            }
        }

        return number_of_zeros;
    }


    /**
     * Determine the number of zero overlaps in the biorthogonal overlap vector.
     *
     * @return The number of zero overlaps.
     */
    template <typename Z = Scalar>
    enable_if_t<std::is_same<Z, complex>::value, int> numberOfZeroOverlaps() const {

        // The count starets at zero.
        int number_of_zeros = 0;

        // Check all overlap values and increase the count if the overlap value is zero.
        for (int i = 0; i < this->biorthogonalOverlaps().rows(); i++) {
            if (std::abs(this->biorthogonalOverlaps()[i].real()) < this->zero_threshold) {
                number_of_zeros += 1;
            }
        }

        return number_of_zeros;
    }


    /**
     * Calculate the reduced overlap. I.e. the biorthogonal overlaps with the zero values get removed and the total remaining overlap is calculated.
     *
     * @return The reduced overlap.
     */
    Scalar reducedOverlap() const {

        // Since the Eigen::vector is essentially a N x 1 matrix, we can remove the rows where the overlap value is zero.
        // The `.removeRows()` method takes a vector and removes all the rows specified by the indices in that vector, which we can calculate with `.zeroOverlapIndices()`.
        auto reduced_overlaps = this->biorthogonalOverlaps();
        reduced_overlaps.removeRows(this->zeroOverlapIndices());

        // We return the product of all remaining overlap elements.
        return reduced_overlaps.prod();
    }


    /**
     * Calculate and return the total overlap value of the biorthogonal coefficients.
     *
     * @return The total overlap value.
     */
    Scalar totalOverlap() const { return this->biorthogonalOverlaps().prod(); }


    /**
     * Determine the indices of the zero overlap values in the biorthogonal overlap vector.
     *
     * @return The indices of the zero overlap values.
     */
    template <typename Z = Scalar>
    enable_if_t<std::is_same<Z, double>::value, std::vector<size_t>> zeroOverlapIndices() const {

        // Initialize the index vector.
        std::vector<size_t> zero_indices {};

        // Check all overlap values and push the index to the index vector if the overlap value is zero.
        for (size_t i = 0; i < this->biorthogonalOverlaps().rows(); i++) {
            if (std::abs(this->biorthogonalOverlaps()[i]) < this->zero_threshold) {
                zero_indices.push_back(i);
            }
        }

        return zero_indices;
    }


    /**
     * Determine the indices of the zero overlap values in the biorthogonal overlap vector.
     *
     * @return The indices of the zero overlap values.
     */
    template <typename Z = Scalar>
    enable_if_t<std::is_same<Z, complex>::value, std::vector<size_t>> zeroOverlapIndices() const {

        // Initialize the index vector.
        std::vector<size_t> zero_indices {};

        // Check all overlap values and push the index to the index vector if the overlap value is zero.
        for (size_t i = 0; i < this->biorthogonalOverlaps().rows(); i++) {
            if (std::abs(this->biorthogonalOverlaps()[i].real()) < this->zero_threshold) {
                zero_indices.push_back(i);
            }
        }

        return zero_indices;
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
     * @note This implementation is based on equation 38b from the 2021 pper by Hugh Burton (https://aip.scitation.org/doi/abs/10.1063/5.0045442).
     */
    DM coDensity(const int k) const {

        // From both the bra and the ket coefficients, we need the k'th column.
        const auto kth_occupied_bra = this->biorthogonalBraExpansion().col(k);
        const auto kth_occupied_ket = this->biorthogonalKetExpansion().col(k);

        // We now perform the product of these two vectors to form the co-density matrix.
        Matrix co_density {this->M, this->M};

        for (int v = 0; v < kth_occupied_ket.rows(); v++) {
            for (int u = 0; u < kth_occupied_bra.rows(); u++) {
                co_density(v, u) = kth_occupied_ket(v) * kth_occupied_bra.conjugate()(u);
            }
        }

        return DM {co_density};
    }


    /**
     * Calculate and return the sum of the zero overlap co-density matrix, the weighted co-density matrix and the zero overlap co-density matrix of the biorthogonalized basis of the bra with itself.
     *
     * @return The sum of the co-density matrices.
     *
     * @note This implementation is based on equation 38d from the 2021 pper by Hugh Burton (https://aip.scitation.org/doi/abs/10.1063/5.0045442).
     */
    DM coDensitySum() const {

        // Initialize the Löwdin pairing basis of the bra with itself.
        Self lp_basis {this->non_orthogonal_state_expansions.first, this->non_orthogonal_state_expansions.first, this->overlap_operator_AO, this->N, 1e-8};

        // Return the sum as described in the function description.
        return lp_basis.zeroOverlapCoDensity() + this->zeroOverlapCoDensity() + this->weightedCoDensity();
    }

    /**
     * Calculate and return the transition one-electron density matrix. The contractions vary depending on the number of zero overlaps (Burton equation 46 and 47).
     * These matrices are then summed and returned.
     *
     * @return The transition one-electron density matrix.
     *
     * @note This implementation is based on equation 45 from the 2021 pper by Hugh Burton (https://aip.scitation.org/doi/abs/10.1063/5.0045442).
     */
    DM transition1DM() const {

        // Initialize a zero matrix.
        Matrix transition_1DM = Matrix::Zero(this->M, this->M);

        // The form of the transition one-electron density matrix depends on the number of zero overlaps.
        if (this->numberOfZeroOverlaps() == 0) {
            transition_1DM += this->reducedOverlap() * this->coDensitySum().matrix();
        } else if (this->numberOfZeroOverlaps() == 1) {
            transition_1DM += this->reducedOverlap() * this->zeroOverlapCoDensity().matrix();
        }

        return DM {transition_1DM.transpose()};
    }


    /**
     * Calculate and return the weighted co-density matrix. It calculates the co-density matrix contributions corresponding to the non-zero overlap values and divides them by that overlap value.
     * These matrices are then summed and returned.
     *
     * @return The weighted co-density matrix.
     *
     * @note This implementation is based on equation 38c from the 2021 paper by Hugh Burton (https://aip.scitation.org/doi/abs/10.1063/5.0045442).
     */
    template <typename Z = Scalar>
    enable_if_t<std::is_same<Z, double>::value, DM> weightedCoDensity() const {


        // Initialize a zero matrix with which we can sum the co-density matrices of the zero overlap indices.
        Matrix weighted_co_density = Matrix::Zero(this->M, this->M);

        // Loop over all zero indices and calculate the zero overlap co-density matrix at each of them. Add them to the total zero overlap co-density matrix.
        for (int i = 0, s = 0; x < this->biorthogonalOverlaps().rows() && s < this->biorthogonalOverlaps().rows(); i++, s++) {
            if (std::abs(this->biorthogonalOverlaps()[s]) > this->zero_threshold) {
                weighted_co_density += (this->coDensity(i).matrix() / this->biorthogonalOverlaps()[s]);
            } else {
                weighted_co_density += this->coDensity(i).matrix();
            }
        }

        return DM {weighted_co_density};
    }


    /**
     * Calculate and return the weighted co-density matrix. It calculates the co-density matrix contributions corresponding to the non-zero overlap values and divides them by that overlap value.
     * These matrices are then summed and returned.
     *
     * @return The weighted co-density matrix.
     *
     * @note This implementation is based on equation 38c from the 2021 paper by Hugh Burton (https://aip.scitation.org/doi/abs/10.1063/5.0045442).
     */
    template <typename Z = Scalar>
    enable_if_t<std::is_same<Z, complex>::value, DM> weightedCoDensity() const {


        // Initialize a zero matrix with which we can sum the co-density matrices of the zero overlap indices.
        Matrix weighted_co_density = Matrix::Zero(this->M, this->M);

        // Loop over all zero indices and calculate the zero overlap co-density matrix at each of them. Add them to the total zero overlap co-density matrix.
        for (int i = 0, s = 0; x < this->biorthogonalOverlaps().rows() && s < this->biorthogonalOverlaps().rows(); i++, s++) {
            if (std::abs(this->biorthogonalOverlaps()[s].real()) > this->zero_threshold) {
                weighted_co_density += (this->coDensity(i).matrix() / this->biorthogonalOverlaps()[s]);
            } else {
                weighted_co_density += this->coDensity(i).matrix();
            }
        }

        return DM {weighted_co_density};
    }


    /**
     * Calculate and return the zero-overlap co-density matrix. It takes the indices of the zero overlaps and calculates the co-density matrix at each of these indices. The sum of those co-density matrices is returned.
     *
     * @return The zero overlap co-density matrix.
     *
     * @note This implementation is based on equation 38a from the 2021 paper by Hugh Burton (https://aip.scitation.org/doi/abs/10.1063/5.0045442).
     */
    DM zeroOverlapCoDensity() const {


        // Initialize a zero matrix with which we can sum the co-density matrices of the zero overlap indices.
        Matrix zero_overlap_co_density = Matrix::Zero(this->M, this->M);

        // Loop over all zero indices and calculate the zero overlap co-density matrix at each of them. Add them to the total zero overlap co-density matrix.
        for (const auto& zero_index : this->zeroOverlapIndices()) {
            zero_overlap_co_density += this->coDensity(zero_index).matrix();
        }

        return DM {zero_overlap_co_density};
    }
};


}  // namespace GQCP
