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


#include "Basis/SpinorBasis/OrbitalSpace.hpp"
#include "Mathematical/Functions/VectorSpaceArithmetic.hpp"
#include "Mathematical/Representation/ImplicitRankFourTensorSlice.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"


namespace GQCP {


/**
 *  The coupled-cluster T2-amplitudes t_{ij}^{ab}.
 * 
 *  In general, these are spinor amplitudes, but they may be used to represent spin-orbital amplitudes as well.
 * 
 *  @param _Scalar          The scalar type that represents one of the amplitudes.
 */
template <typename _Scalar>
class T2Amplitudes:
    public VectorSpaceArithmetic<T2Amplitudes<_Scalar>, _Scalar> {

public:
    // The scalar type that represents one of the amplitudes.
    using Scalar = _Scalar;

    // The type of 'this'.
    using Self = T2Amplitudes<Scalar>;

private:
    // The orbital space which encapsulates the occupied-virtual separation.
    OrbitalSpace orbital_space;

    // The T2-amplitudes as an implicit tensor, implementing intuitive operator(i,j,a,b) calls.
    ImplicitRankFourTensorSlice<Scalar> t;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  Construct T2-amplitudes given their representation as an `ImplicitRankFourTensorSlice` and explicit occupied-virtual orbital space.
     * 
     *  @param t                            The T2-amplitudes as an implicit tensor, implementing intuitive operator(i,j,a,b) calls.
     *  @param orbital_space                The orbital space which encapsulates the occupied-virtual separation.
     */
    T2Amplitudes(const ImplicitRankFourTensorSlice<Scalar>& t, const OrbitalSpace& orbital_space) :
        orbital_space {orbital_space},
        t {t} {}


    /**
     *  Construct the T2-amplitudes given their representation as an `ImplicitRankFourTensorSlice` and an implicit occupied-virtual orbital space determined by the given number of occupied orbitals and total number of orbitals.
     * 
     *  @param t                The T2-amplitudes as an implicit tensor, implementing intuitive operator(i,j,a,b) calls.
     *  @param N                The number of occupied orbitals.
     *  @param M                The total number of orbitals.
     */
    T2Amplitudes(const ImplicitRankFourTensorSlice<Scalar>& t, const size_t N, const size_t M) :
        T2Amplitudes(t, OrbitalSpace::Implicit({{OccupationType::k_occupied, N}, {OccupationType::k_virtual, M - N}})) {}


    /*
     *  MARK: Named constructors
     */

    /**
     *  Create perturbative T2-amplitudes.
     * 
     *  @param f                            The (inactive) Fock matrix.
     *  @param V_A                          The antisymmetrized two-electron integrals (in physicist's notation).
     *  @param orbital_space                The orbital space which encapsulates the occupied-virtual separation.
     * 
     *  @return T2-amplitudes calculated from an initial perturbative result.
     */
    static T2Amplitudes<Scalar> Perturbative(const SquareMatrix<Scalar>& f, const SquareRankFourTensor<Scalar>& V_A, const OrbitalSpace& orbital_space) {

        // Zero-initialize a tensor representation for the (occupied-occupied-virtual-virtual) T2-amplitudes t_{ij}^{ab}.
        auto t2 = orbital_space.initializeRepresentableObjectFor<Scalar>(OccupationType::k_occupied, OccupationType::k_occupied, OccupationType::k_virtual, OccupationType::k_virtual);


        // Provide the perturbative T2-amplitudes, by setting the RHS amplitudes in Stanton1991, equation (2) to zero.
        for (const auto& i : orbital_space.indices(OccupationType::k_occupied)) {
            for (const auto& j : orbital_space.indices(OccupationType::k_occupied)) {
                for (const auto& a : orbital_space.indices(OccupationType::k_virtual)) {
                    for (const auto& b : orbital_space.indices(OccupationType::k_virtual)) {
                        const auto denominator = f(i, i) + f(j, j) - f(a, a) - f(b, b);

                        t2(i, j, a, b) = V_A(i, j, a, b) / denominator;
                    }
                }
            }
        }

        return T2Amplitudes<Scalar>(t2, orbital_space);
    }


    /*
     *  MARK: Access
     */

    /**
     *  Access one of the T2-amplitudes.
     * 
     *  @param i            An occupied index.
     *  @param j            An occupied index.
     *  @param a            A virtual index.
     *  @param b            A virtual index.
     * 
     *  @return A read-only reference to the T2-amplitude t_{ij}^{ab}.
     */
    Scalar operator()(const size_t i, const size_t j, const size_t a, const size_t b) const { return this->t(i, j, a, b); }

    /**
     *  Access one of the T2-amplitudes.
     * 
     *  @param i            An occupied index.
     *  @param j            An occupied index.
     *  @param a            A virtual index.
     *  @param b            A virtual index.
     * 
     *  @return A writable reference to the T2-amplitude corresponding to t_{ij}^{ab}
     */
    Scalar& operator()(const size_t i, const size_t j, const size_t a, const size_t b) { return this->t(i, j, a, b); }

    /**
     *  @return The T2-amplitudes as an `ImplicitRankFourTensorSlice`.
     */
    const ImplicitRankFourTensorSlice<Scalar>& asImplicitRankFourTensorSlice() const { return this->t; }

    /**
     *  @return The orbital space for these T2-amplitudes, which encapsulates the occupied-virtual separation.
     */
    const OrbitalSpace& orbitalSpace() const { return this->orbital_space; }


    /*
     *  MARK: Linear algebra
     */

    /**
     *  @return The Frobenius norm of these T2-amplitudes.
     */
    Scalar norm() const { return this->asImplicitRankFourTensorSlice().asMatrix().norm(); }


    /*
     *  MARK: Vector space arithmetic
     */

    /**
     *  Addition-assignment.
     */
    Self& operator+=(const Self& rhs) override {

        // Prepare some variables.
        const auto& index_maps = this->asImplicitRankFourTensorSlice().indexMaps();

        // Add the tensor representations.
        const Tensor<Scalar, 4> t_sum_dense = this->asImplicitRankFourTensorSlice().asTensor().Eigen() + rhs.asImplicitRankFourTensorSlice().asTensor().Eigen();
        const ImplicitRankFourTensorSlice<Scalar> t_sum_slice {index_maps, t_sum_dense};

        this->t = t_sum_slice;

        return *this;
    }


    /**
     *  Scalar multiplication-assignment.
     */
    Self& operator*=(const Scalar& a) override {

        // Prepare some variables.
        const auto& index_maps = this->asImplicitRankFourTensorSlice().indexMaps();

        // Multiply the tensor representation.
        const Tensor<Scalar, 4> t_multiplied_dense = a * this->asImplicitRankFourTensorSlice().asTensor().Eigen();
        const ImplicitRankFourTensorSlice<Scalar> t_multiplied_slice {index_maps, t_multiplied_dense};

        this->t = t_multiplied_slice;

        return *this;
    }
};

}  // namespace GQCP
