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
#include "Mathematical/Representation/ImplicitMatrixSlice.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "Utilities/type_traits.hpp"


namespace GQCP {


/**
 *  The coupled-cluster T1 amplitudes t_i^a.
 * 
 *  In general, these are spinor amplitudes, but they may be used to represent spin-orbital amplitudes as well.
 * 
 *  @param _Scalar          The scalar type that represents one of the amplitudes.
 */
template <typename _Scalar>
class T1Amplitudes:
    public VectorSpaceArithmetic<T1Amplitudes<_Scalar>, _Scalar> {

public:
    // The scalar type that represents one of the amplitudes.
    using Scalar = _Scalar;

    // The type of 'this'.
    using Self = T1Amplitudes<Scalar>;


private:
    // The orbital space which encapsulates the occupied-virtual separation.
    OrbitalSpace orbital_space;

    // The T1-amplitudes as an implicit matrix, implementing intuitive operator(i,a) calls.
    ImplicitMatrixSlice<Scalar> t;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  Construct T1-amplitudes given their representation as a ImplicitMatrixSlice and explicit occupied-virtual orbital space.
     * 
     *  @param t                            The T1-amplitudes as an implicit matrix, implementing intuitive operator(i,a) calls.
     *  @param orbital_space                The orbital space which encapsulates the occupied-virtual separation.
     */
    T1Amplitudes(const ImplicitMatrixSlice<Scalar>& t, const OrbitalSpace& orbital_space) :
        orbital_space {orbital_space},
        t {t} {}


    /**
     *  Construct the T1-amplitudes given their representation as an `ImplicitMatrixSlice` and an implicit occupied-virtual orbital space determined by the given number of occupied orbitals and total number of orbitals.
     * 
     *  @param t                The T1-amplitudes as an implicit matrix, implementing intuitive operator(i,a) calls.
     *  @param N                The number of occupied orbitals.
     *  @param M                The total number of orbitals.
     */
    T1Amplitudes(const ImplicitMatrixSlice<Scalar>& t, const size_t N, const size_t M) :
        T1Amplitudes(t, OrbitalSpace::Implicit({{OccupationType::k_occupied, N}, {OccupationType::k_virtual, M - N}})) {}


    /*
     *  MARK: Named constructors
     */

    /**
     *  Create perturbative T1-amplitudes.
     * 
     *  @param f                            The (inactive) Fock matrix.
     *  @param orbital_space                The orbital space which covers the occupied-virtual separation.
     * 
     *  @return T1-amplitudes calculated from an initial perturbative result.
     */
    static T1Amplitudes<Scalar> Perturbative(const SquareMatrix<Scalar>& f, const OrbitalSpace& orbital_space) {

        // Zero-initialize a matrix representation for the (occupied-virtual) T1-amplitudes t_i^a.
        auto t1 = orbital_space.initializeRepresentableObjectFor<Scalar>(OccupationType::k_occupied, OccupationType::k_virtual);

        // Provide the perturbative T1-amplitudes, by setting the RHS amplitudes in Stanton1991, equation (1) to zero.
        for (const auto& i : orbital_space.indices(OccupationType::k_occupied)) {
            for (const auto& a : orbital_space.indices(OccupationType::k_virtual)) {
                const auto denominator = f(i, i) - f(a, a);

                t1(i, a) = f(i, a) / denominator;
            }
        }

        return T1Amplitudes<Scalar>(t1, orbital_space);
    }


    /*
     *  MARK: Access.
     */

    /**
     *  Access one of the T1-amplitudes.
     * 
     *  @param i            An occupied index.
     *  @param a            A virtual index.
     * 
     *  @return A read-only reference to the T1-amplitude t_i^a.
     */
    Scalar operator()(const size_t i, const size_t a) const { return this->t(i, a); }

    /**
     *  Access one of the T1-amplitudes.
     * 
     *  @param i            An occupied index.
     *  @param a            A virtual index.
     * 
     *  @return A writable reference to the T1-amplitude t_i^a.
     */
    Scalar& operator()(const size_t i, const size_t a) { return this->t(i, a); }

    /**
     *  @return The T1-amplitudes as an `ImplicitMatrixSlice`.
     */
    const ImplicitMatrixSlice<Scalar>& asImplicitMatrixSlice() const { return this->t; }

    /**
     *  @return The orbital space for these T1-amplitudes, which encapsulates the occupied-virtual separation.
     */
    const OrbitalSpace& orbitalSpace() const { return this->orbital_space; }


    /*
     *  MARK: Linear algebra
     */

    /**
     *  @return The Frobenius norm of these amplitudes. 
     */
    Scalar norm() const { return this->asImplicitMatrixSlice().asMatrix().norm(); }


    /*
     *  MARK: Vector space arithmetic
     */

    /**
     *  Addition-assignment.
     */
    Self& operator+=(const Self& rhs) override {

        // Prepare some variables.
        const auto& row_map = this->asImplicitMatrixSlice().rowIndexMap();
        const auto& col_map = this->asImplicitMatrixSlice().columnIndexMap();

        // Sum the matrix representations.
        const MatrixX<Scalar> t_sum_dense = this->asImplicitMatrixSlice().asMatrix() + rhs.asImplicitMatrixSlice().asMatrix();
        const ImplicitMatrixSlice<Scalar> t_sum_slice {row_map, col_map, t_sum_dense};

        this->t = t_sum_slice;

        return *this;
    }


    /**
     *  Scalar multiplication-assignment.
     */
    Self& operator*=(const Scalar& a) override {

        // Prepare some variables.
        const auto& row_map = this->asImplicitMatrixSlice().rowIndexMap();
        const auto& col_map = this->asImplicitMatrixSlice().columnIndexMap();

        // Multiply the matrix representation.
        const MatrixX<Scalar> t_multiplied_dense = a * this->asImplicitMatrixSlice().asMatrix();
        const ImplicitMatrixSlice<Scalar> t_multiplied_slice {row_map, col_map, t_multiplied_dense};

        this->t = t_multiplied_slice;

        return *this;
    }
};

}  // namespace GQCP
