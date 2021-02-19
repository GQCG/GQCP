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
 *  A base class for the representation of doubles amplitudes. It is mainly used as a class that implements common behavior for the classes `T2Amplitudes` and `L2Amplitudes`.
 * 
 *  @param _Scalar                  The scalar type of one of the amplitudes.
 *  @param _DerivedAmplitudes       The type of amplitudes class that derives from this class, enabling CRTP and compile-time polymorphism.
 */
template <typename _Scalar, typename _DerivedAmplitudes>
class DoublesAmplitudes:
    public VectorSpaceArithmetic<_DerivedAmplitudes, _Scalar> {

public:
    // The scalar type of one of the amplitudes.
    using Scalar = _Scalar;

    // The type of amplitudes class that derives from this class, enabling CRTP and compile-time polymorphism.
    using DerivedAmplitudes = _DerivedAmplitudes;


private:
    // The orbital space which encapsulates the occupied-virtual separation.
    OrbitalSpace orbital_space;

    // The doubles amplitudes as an implicit tensor, implementing intuitive operator(i,j,a,b) calls
    ImplicitRankFourTensorSlice<Scalar> v;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  Construct doubles amplitudes given their representation as an `ImplicitRankFourTensorSlice` and an explicit occupied-virtual orbital space.
     * 
     *  @param v                        The doubles amplitudes as an implicit matrix, implementing intuitive operator(i,a) calls.
     *  @param orbital_space            The orbital space which encapsulates the occupied-virtual separation.
     */
    DoublesAmplitudes(const ImplicitRankFourTensorSlice<Scalar>& v, const OrbitalSpace& orbital_space) :
        orbital_space {orbital_space},
        v {v} {}


    /**
     *  Construct doubles amplitudes given their representation as an `ImplicitRankFourTensorSlice` and an implicit occupied-virtual orbital space determined by the given number of occupied orbitals and total number of orbitals.
     * 
     *  @param v                The doubles amplitudes as an implicit matrix, implementing intuitive operator(i,a) calls.
     *  @param N                The number of occupied orbitals.
     *  @param M                The total number of orbitals.
     */
    DoublesAmplitudes(const ImplicitRankFourTensorSlice<Scalar>& v, const size_t N, const size_t M) :
        DoublesAmplitudes(v, OrbitalSpace::Implicit({{OccupationType::k_occupied, N}, {OccupationType::k_virtual, M - N}})) {}


    /*
     *  MARK: Access
     */

    /**
     *  Access one of the doubles amplitudes.
     * 
     *  @param i            An occupied index.
     *  @param j            An occupied index.
     *  @param a            A virtual index.
     *  @param b            A virtual index.
     * 
     *  @return A read-only reference to the doubles amplitude u_{ij}^{ab}.
     */
    Scalar operator()(const size_t i, const size_t j, const size_t a, const size_t b) const { return this->v(i, j, a, b); }

    /**
     *  Access one of the doubles amplitudes.
     * 
     *  @param i            An occupied index.
     *  @param j            An occupied index.
     *  @param a            A virtual index.
     *  @param b            A virtual index.
     * 
     *  @return A writable reference to the doubles amplitude u_{ij}^{ab}.
     */
    Scalar& operator()(const size_t i, const size_t j, const size_t a, const size_t b) { return this->v(i, j, a, b); }

    /**
     *  @return These doubles amplitudes as an `ImplicitRankFourTensorSlice`.
     */
    const ImplicitRankFourTensorSlice<Scalar>& asImplicitRankFourTensorSlice() const { return this->v; }

    /**
     *  @return The orbital space which encapsulates the occupied-virtual separation.
     */
    const OrbitalSpace& orbitalSpace() const { return this->orbital_space; }


    /*
     *  MARK: Linear algebra
     */

    /**
     *  @return The Frobenius norm of these doubles amplitudes.
     */
    Scalar norm() const { return this->asImplicitRankFourTensorSlice().asMatrix().norm(); }


    /*
     *  MARK: Vector space arithmetic
     */

    /**
     *  Addition-assignment.
     */
    DerivedAmplitudes& operator+=(const DerivedAmplitudes& rhs) override {

        // Prepare some variables.
        const auto& index_maps = this->asImplicitRankFourTensorSlice().indexMaps();

        // Add the tensor representation.
        const Tensor<Scalar, 4> v_sum_dense = this->asImplicitRankFourTensorSlice().asTensor().Eigen() + rhs.asImplicitRankFourTensorSlice().asTensor().Eigen();
        const ImplicitRankFourTensorSlice<Scalar> v_sum_slice {index_maps, v_sum_dense};

        this->v = v_sum_slice;

        return static_cast<DerivedAmplitudes&>(*this);
    }


    /**
     *  Scalar multiplication-assignment.
     */
    DerivedAmplitudes& operator*=(const Scalar& a) override {

        // Prepare some variables.
        const auto& index_maps = this->asImplicitRankFourTensorSlice().indexMaps();

        // Multiply the tensor representation.
        const Tensor<Scalar, 4> v_multiplied_dense = a * this->asImplicitRankFourTensorSlice().asTensor().Eigen();
        const ImplicitRankFourTensorSlice<Scalar> v_multiplied_slice {index_maps, v_multiplied_dense};

        this->v = v_multiplied_slice;

        return static_cast<DerivedAmplitudes&>(*this);
    }
};

}  // namespace GQCP
