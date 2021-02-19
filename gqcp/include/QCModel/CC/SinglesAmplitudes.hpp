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
#include "Utilities/type_traits.hpp"


namespace GQCP {


/**
 *  A base class for the representation of singles amplitudes. It is mainly used as a class that implements common behavior for the classes `T1Amplitudes` and `L1Amplitudes`.
 * 
 *  @param _Scalar                  The scalar type of one of the amplitudes.
 *  @param _DerivedAmplitudes       The type of amplitudes class that derives from this class, enabling CRTP and compile-time polymorphism.
 */
template <typename _Scalar, typename _DerivedAmplitudes>
class SinglesAmplitudes:
    public VectorSpaceArithmetic<_DerivedAmplitudes, _Scalar> {

public:
    // The scalar type of one of the amplitudes.
    using Scalar = _Scalar;

    // The type of amplitudes class that derives from this class, enabling CRTP and compile-time polymorphism.
    using DerivedAmplitudes = _DerivedAmplitudes;


private:
    // The orbital space which encapsulates the occupied-virtual separation.
    OrbitalSpace orbital_space;

    // The singles amplitudes as an implicit matrix, implementing intuitive operator(i,a) calls.
    ImplicitMatrixSlice<Scalar> u;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  Construct singles amplitudes given their representation as an `ImplicitMatrixSlice` and an explicit occupied-virtual orbital space.
     * 
     *  @param u                        The singles amplitudes as an implicit matrix, implementing intuitive operator(i,a) calls.
     *  @param orbital_space            The orbital space which encapsulates the occupied-virtual separation.
     */
    SinglesAmplitudes(const ImplicitMatrixSlice<Scalar>& u, const OrbitalSpace& orbital_space) :
        orbital_space {orbital_space},
        u {u} {}


    /**
     *  Construct singles amplitudes given their representation as an `ImplicitMatrixSlice` and an implicit occupied-virtual orbital space determined by the given number of occupied orbitals and total number of orbitals.
     * 
     *  @param u                The singles amplitudes as an implicit matrix, implementing intuitive operator(i,a) calls.
     *  @param N                The number of occupied orbitals.
     *  @param M                The total number of orbitals.
     */
    SinglesAmplitudes(const ImplicitMatrixSlice<Scalar>& u, const size_t N, const size_t M) :
        SinglesAmplitudes(u, OrbitalSpace::Implicit({{OccupationType::k_occupied, N}, {OccupationType::k_virtual, M - N}})) {}


    /*
     *  MARK: Access
     */

    /**
     *  Access one of the singles amplitudes.
     * 
     *  @param p            The subscript index.
     *  @param q            The superscript index.
     * 
     *  @return A read-only reference to the singles amplitude u_p^q.
     */
    Scalar operator()(const size_t p, const size_t q) const { return this->u(p, q); }

    /**
     *  Access one of the singles amplitudes.
     * 
     *  @param p            The subscript index.
     *  @param q            The superscript index.
     * 
     *  @return A writable reference to the singles amplitude u_p^q.
     */
    Scalar& operator()(const size_t p, const size_t q) { return this->u(p, q); }

    /**
     *  @return These singles amplitudes as an `ImplicitMatrixSlice`.
     */
    const ImplicitMatrixSlice<Scalar>& asImplicitMatrixSlice() const { return this->u; }

    /**
     *  @return The orbital space which encapsulates the occupied-virtual separation.
     */
    const OrbitalSpace& orbitalSpace() const { return this->orbital_space; }


    /*
     *  MARK: Linear algebra
     */

    /**
     *  @return The Frobenius norm of these singles amplitudes.
     */
    Scalar norm() const { return this->asImplicitMatrixSlice().asMatrix().norm(); }


    /*
     *  MARK: Vector space arithmetic
     */

    /**
     *  Addition-assignment.
     */
    DerivedAmplitudes& operator+=(const DerivedAmplitudes& rhs) override {

        // Prepare some variables.
        const auto& row_map = this->asImplicitMatrixSlice().rowIndexMap();
        const auto& col_map = this->asImplicitMatrixSlice().columnIndexMap();

        // Add the matrix representation.
        const MatrixX<Scalar> u_sum_dense = this->asImplicitMatrixSlice().asMatrix() + rhs.asImplicitMatrixSlice().asMatrix();
        const ImplicitMatrixSlice<Scalar> u_sum_slice {row_map, col_map, u_sum_dense};

        this->u = u_sum_slice;

        return static_cast<DerivedAmplitudes&>(*this);
    }


    /**
     *  Scalar multiplication-assignment.
     */
    DerivedAmplitudes& operator*=(const Scalar& a) override {

        // Prepare some variables.
        const auto& row_map = this->asImplicitMatrixSlice().rowIndexMap();
        const auto& col_map = this->asImplicitMatrixSlice().columnIndexMap();

        // Multiply the matrix representation.
        const MatrixX<Scalar> u_multiplied_dense = a * this->asImplicitMatrixSlice().asMatrix();
        const ImplicitMatrixSlice<Scalar> u_multiplied_slice {row_map, col_map, u_multiplied_dense};

        this->u = u_multiplied_slice;

        return static_cast<DerivedAmplitudes&>(*this);
    }
};

}  // namespace GQCP
