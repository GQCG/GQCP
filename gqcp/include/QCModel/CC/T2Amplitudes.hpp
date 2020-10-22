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
 *  The coupled-cluster T2-amplitudes t_{ij}^{ab}. According to context, this class may either represent restricted (i.e. spatial-orbital) amplitudes, or generalized (spinor) amplitudes.
 */
template <typename _Scalar>
class T2Amplitudes:
    public VectorSpaceArithmetic<T2Amplitudes<_Scalar>, _Scalar> {

public:
    using Scalar = _Scalar;
    using Self = T2Amplitudes<Scalar>;

private:
    OrbitalSpace orbital_space;  // the orbital space which covers the occupied-virtual separation

    ImplicitRankFourTensorSlice<Scalar> t;  // the T2-amplitudes as a block rank-four tensor, implementing easy operator(i,j,a,b) calls


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  Construct T2-amplitudes given their representation as a ImplicitRankFourTensorSlice and and explicit occupied-virtual orbital space.
     * 
     *  @param t                            the T2-amplitudes as a block matrix, implementing easy operator(i,j,a,b) calls
     *  @param orbital_space                the orbital space which covers the occupied-virtual separation, indicating which indices in the block matrix are occupied and which are virtual
     */
    T2Amplitudes(const ImplicitRankFourTensorSlice<Scalar>& t, const OrbitalSpace& orbital_space) :
        orbital_space {orbital_space},
        t {t} {}


    /**
     *  Construct the T2-amplitudes given their representation as a ImplicitRankFourTensorSlice and an implicit occupied-virtual orbital space determined by the given number of occupied orbitals and total number of orbitals.
     * 
     *  @param t                the T2-amplitudes as a block matrix, implementing easy operator(i,j,a,b) calls
     *  @param N                the number of occupied orbitals
     *  @param M                the total number of orbitals
     */
    T2Amplitudes(const ImplicitRankFourTensorSlice<Scalar>& t, const size_t N, const size_t M) :
        T2Amplitudes(t, OrbitalSpace::Implicit({{OccupationType::k_occupied, N}, {OccupationType::k_virtual, M - N}})) {}


    /*
     *  NAMED CONSTRUCTORS
     */

    /**
     *  Create perturbative T2-amplitudes using an explicit orbital space.
     * 
     *  @param f                            the (inactive) Fock matrix
     *  @param V_A                          the antisymmetrized two-electron integrals (in physicist's notation)
     *  @param orbital_space                the orbital space which covers the occupied-virtual separation
     * 
     *  @return T2-amplitudes calculated from an initial perturbative result
     */
    static T2Amplitudes<Scalar> Perturbative(const SquareMatrix<Scalar>& f, const QCRankFourTensor<Scalar>& V_A, const OrbitalSpace& orbital_space) {

        // Zero-initialize a tensor representation for the (occupied-occupied-virtual-virtual) T2-amplitudes t_{ij}^{ab}.
        auto t2 = orbital_space.initializeRepresentableObjectFor<double>(OccupationType::k_occupied, OccupationType::k_occupied, OccupationType::k_virtual, OccupationType::k_virtual);


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


    /**
     *  Create perturbative T2-amplitudes using an implicit occupied-virtual orbital space determined by the given number of occupied orbitals and total number of orbitals.
     * 
     *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal spinor basis
     *  @param N                            the number of occupied orbitals
     *  @param M                            the total number of orbitals
     * 
     *  @return T2-amplitudes calculated from an initial perturbative result
     */
    static T2Amplitudes<Scalar> Perturbative(const RSQHamiltonian<Scalar>& sq_hamiltonian, const size_t N, const size_t M) {

        // Create the implicit orbital space for N occupied orbitals and M total orbitals.
        const auto orbital_space = OrbitalSpace::Implicit({{OccupationType::k_occupied, N}, {OccupationType::k_virtual, M - N}});

        return T2Amplitudes<Scalar>::Perturbative(sq_hamiltonian, orbital_space);
    }


    /*
     *  OPERATORS
     */

    /**
     *  Access one of the T2-amplitudes.
     * 
     *  @param i            an occupied index
     *  @param j            an occupied index
     *  @param a            a virtual index
     *  @param b            a virtual index
     * 
     *  @return the read-only T2-amplitude corresponding to t_{ij}^{ab}
     */
    Scalar operator()(const size_t i, const size_t j, const size_t a, const size_t b) const { return this->t(i, j, a, b); }

    /**
     *  Access one of the T2-amplitudes.
     * 
     *  @param i            an occupied index
     *  @param j            an occupied index
     *  @param a            a virtual index
     *  @param b            a virtual index
     * 
     *  @return the writable T2-amplitude corresponding to t_{ij}^{ab}
     */
    Scalar& operator()(const size_t i, const size_t j, const size_t a, const size_t b) { return this->t(i, j, a, b); }

    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return the T2-amplitudes as a ImplicitRankFourTensorSlice
     */
    const ImplicitRankFourTensorSlice<Scalar>& asImplicitRankFourTensorSlice() const { return this->t; }

    /**
     *  @return the Frobenius norm of these T2-amplitudes
     */
    Scalar norm() const { return this->asImplicitRankFourTensorSlice().asMatrix().norm(); }

    /**
     *  @return the orbital space for these T2-amplitudes, which covers the occupied-virtual separation
     */
    const OrbitalSpace& orbitalSpace() const { return this->orbital_space; }


    /*
     *  MARK: Vector space arithmetic
     */

    /**
     *  Addition-assignment.
     */
    Self& operator+=(const Self& rhs) override {

        // Prepare some variables.
        const auto& index_maps = this->asImplicitRankFourTensorSlice().indexMaps();

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

        const Tensor<Scalar, 4> t_multiplied_dense = a * this->asImplicitRankFourTensorSlice().asTensor().Eigen();
        const ImplicitRankFourTensorSlice<Scalar> t_multiplied_slice {index_maps, t_multiplied_dense};

        this->t = t_multiplied_slice;

        return *this;
    }
};

}  // namespace GQCP
