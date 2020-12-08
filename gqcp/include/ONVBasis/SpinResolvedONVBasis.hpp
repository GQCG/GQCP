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


#include "ONVBasis/SpinUnresolvedONVBasis.hpp"
#include "Operator/SecondQuantized/MixedUSQTwoElectronOperatorComponent.hpp"
#include "Operator/SecondQuantized/ModelHamiltonian/HubbardHamiltonian.hpp"
#include "Operator/SecondQuantized/RSQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/RSQTwoElectronOperator.hpp"
#include "QuantumChemical/SpinResolvedBase.hpp"

#include <functional>


namespace GQCP {


/**
 *  A full spin-resolved ONV basis.
 */
class SpinResolvedONVBasis:
    public SpinResolvedBase<SpinUnresolvedONVBasis, SpinResolvedONVBasis> {
public:
    // The type component this spin resolved object is made of.
    using ComponentType = typename SpinResolvedBase<SpinUnresolvedONVBasis, SpinResolvedONVBasis>::Of;

private:
    // A vector of sparse matrices containing the one-electron coupling elements for the alpha ONV basis. See also `calculateOneElectronCouplings`.
    std::vector<Eigen::SparseMatrix<double>> alpha_couplings;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  @param K            The number of alpha or beta spin-orbitals.
     *  @param N_alpha      The number of alpha electrons, i.e. the number of occupied alpha spin-orbitals.
     *  @param N_beta       The number of beta electrons, i.e. the number of occupied beta spin-orbitals.
     */
    SpinResolvedONVBasis(const size_t K, const size_t N_alpha, const size_t N_beta);


    /*
     *  MARK: General information
     */

    /**
     *  Calculate the dimension of a spin-resolved ONV basis with a given number of orbitals and electrons.
     * 
     *  @param K            The number of alpha or beta spin-orbitals.
     *  @param N_alpha      The number of alpha electrons, i.e. the number of occupied alpha spin-orbitals.
     *  @param N_beta       The number of beta electrons, i.e. the number of occupied beta spin-orbitals.
     *
     *  @return The dimension of a spin-resolved ONV basis.
     */
    static size_t calculateDimension(const size_t K, const size_t N_alpha, const size_t N_beta);

    /**
     *  @return The dimension of this ONV basis.
     */
    size_t dimension() const;

    /**
     *  @return The number of alpha or beta spin-orbitals.
     */
    size_t numberOfOrbitals() const { return this->alpha().numberOfOrbitals(); }


    /*
     *  MARK: Couplings
     */

    /**
     *  @return A vector of sparse matrices containing the one-electron coupling elements for the alpha ONV basis. See also `calculateOneElectronCouplings`.
     */
    const std::vector<Eigen::SparseMatrix<double>>& alphaCouplings() const { return alpha_couplings; }

    /**
     *  Calculate the one-electron operator intermediate that is required for the calculation of "theta(pq)" in Helgaker, Jørgensen, Olsen (2000). It is a partitioning of the mixed component of the unrestricted two-electron operator g(ab)_{pqrs}, resulting in a one-electron operator t(b)_{rs}.
     *
     *  @param p            The first index of the two-electron operator.
     *  @param q            The second index of the two-electron operator.
     *  @param g_ab_op      The two-electron operator.
     *
     *  @return The intermediate one-electron operator that is required for the calculation of "theta(pq)".
     */
    ScalarUSQOneElectronOperatorComponent<double> calculateOneElectronPartition(const size_t p, const size_t q, const ScalarMixedUSQTwoElectronOperatorComponent<double>& g_ab_op) const;


    /*
     *  MARK: Address calculations
     */

    /**
     *  Calculate the compound address of an ONV represented by the two given alpha- and beta-addresses.
     * 
     *  @param I_alpha              The alpha-address.
     *  @param I_beta               The beta-address.
     * 
     *  @return The compound address of an ONV represented by the two given alpha- and beta-addresses.
     */
    size_t compoundAddress(const size_t I_alpha, const size_t I_beta) const;


    /*
     *  MARK: Iterations
     */

    /**
     *  Iterate over all ONVs (implicitly, by resolving in their spin components) in this ONV basis and apply the given callback function.
     * 
     *  @param callback             The function to be applied in every iteration. Its arguments are two pairs of spin-unresolved ONVs and their corresponding addresses, where the first two arguments are related to alpha-spin. The last two arguments are related to beta-spin.
     */
    void forEach(const std::function<void(const SpinUnresolvedONV&, const size_t, const SpinUnresolvedONV&, const size_t)>& callback) const;


    /*
     *  MARK: Dense restricted operator evaluations
     */

    /**
     *  Calculate the dense matrix representation of a restricted one-electron operator in this ONV basis.
     *
     *  @param f                A restricted one-electron operator expressed in an orthonormal orbital basis.
     *
     *  @return A dense matrix represention of the one-electron operator.
     */
    SquareMatrix<double> evaluateOperatorDense(const ScalarRSQOneElectronOperator<double>& f) const;

    /**
     *  Calculate the dense matrix representation of a restricted two-electron operator in this ONV basis.
     *
     *  @param g                A restricted two-electron operator expressed in an orthonormal orbital basis.
     *
     *  @return A dense matrix represention of the two-electron operator.
     */
    SquareMatrix<double> evaluateOperatorDense(const ScalarRSQTwoElectronOperator<double>& g) const;

    /**
     *  Calculate the dense matrix representation of a restricted Hamiltonian in this ONV basis.
     *
     *  @param hamiltonian      A restricted Hamiltonian expressed in an orthonormal orbital basis.
     *
     *  @return A dense matrix represention of the Hamiltonian.
     */
    SquareMatrix<double> evaluateOperatorDense(const RSQHamiltonian<double>& hamiltonian) const;

    /**
     *  Calculate the dense matrix representation of a Hubbard Hamiltonian in this ONV basis.
     *
     *  @param hamiltonian      A Hubbard Hamiltonian expressed in an orthonormal orbital basis.
     *
     *  @return A dense matrix represention of the Hamiltonian.
     */
    SquareMatrix<double> evaluateOperatorDense(const HubbardHamiltonian<double>& hamiltonian) const;


    /*
     *  MARK: Diagonal restricted operator evaluations
     */

    /**
     *  Calculate the diagonal of the matrix representation of a restricted one-electron operator in this ONV basis.
     *
     *  @param f_op             A restricted one-electron operator expressed in an orthonormal orbital basis.
     *
     *  @return The diagonal of the dense matrix represention of the one-electron operator.
     */
    VectorX<double> evaluateOperatorDiagonal(const ScalarRSQOneElectronOperator<double>& f_op) const;

    /**
     *  Calculate the diagonal of the matrix representation of a restricted two-electron operator in this ONV basis.
     *
     *  @param g                A restricted two-electron operator expressed in an orthonormal orbital basis.
     *
     *  @return The diagonal of the dense matrix represention of the two-electron operator.
     */
    VectorX<double> evaluateOperatorDiagonal(const ScalarRSQTwoElectronOperator<double>& g) const;

    /**
     *  Calculate the diagonal of the dense matrix representation of a restricted Hamiltonian in this ONV basis.
     *
     *  @param hamiltonian      A restricted Hamiltonian expressed in an orthonormal orbital basis.
     *
     *  @return The diagonal of the dense matrix represention of the Hamiltonian.
     */
    VectorX<double> evaluateOperatorDiagonal(const RSQHamiltonian<double>& hamiltonian) const;

    /**
     *  Calculate the diagonal of the dense matrix representation of a Hubbard Hamiltonian in this ONV basis.
     *
     *  @param hamiltonian      A Hubbard Hamiltonian expressed in an orthonormal orbital basis.
     *
     *  @return The diagonal of the dense matrix represention of the Hamiltonian.
     */
    VectorX<double> evaluateOperatorDiagonal(const HubbardHamiltonian<double>& hamiltonian) const;


    /*
     *  MARK: Restricted matrix-vector product evaluations
     */

    /**
     *  Calculate the matrix-vector product of (the matrix representation of) a restricted one-electron operator with the given coefficient vector.
     *
     *  @param f                A restricted one-electron operator expressed in an orthonormal orbital basis.
     *  @param x                The coefficient vector of a linear expansion.
     *
     *  @return The coefficient vector of the linear expansion after being acted on with the given (matrix representation of) the one-electron operator.
     */
    VectorX<double> evaluateOperatorMatrixVectorProduct(const ScalarRSQOneElectronOperator<double>& f, const VectorX<double>& x) const;

    /**
     *  Calculate the matrix-vector product of (the matrix representation of) a restricted two-electron operator with the given coefficient vector.
     *
     *  @param g                A restricted two-electron operator expressed in an orthonormal orbital basis.
     *  @param x                The coefficient vector of a linear expansion.
     *
     *  @return The coefficient vector of the linear expansion after being acted on with the given (matrix representation of) the two-electron operator.
     */
    VectorX<double> evaluateOperatorMatrixVectorProduct(const ScalarRSQTwoElectronOperator<double>& g, const VectorX<double>& x) const;

    /**
     *  Calculate the matrix-vector product of (the matrix representation of) a restricted Hamiltonian with the given coefficient vector.
     *
     *  @param hamiltonian      A restricted Hamiltonian expressed in an orthonormal orbital basis.
     *  @param x                The coefficient vector of a linear expansion.
     *
     *  @return The coefficient vector of the linear expansion after being acted on with the given (matrix representation of) the Hamiltonian.
     */
    VectorX<double> evaluateOperatorMatrixVectorProduct(const RSQHamiltonian<double>& hamiltonian, const VectorX<double>& x) const;

    /**
     *  Calculate the matrix-vector product of (the matrix representation of) a Hubbard Hamiltonian with the given coefficient vector.
     *
     *  @param hamiltonian      A Hubbard Hamiltonian expressed in an orthonormal orbital basis.
     *  @param x                The coefficient vector of a linear expansion.
     *
     *  @return The coefficient vector of the linear expansion after being acted on with the given (matrix representation of) the Hamiltonian.
     */
    VectorX<double> evaluateOperatorMatrixVectorProduct(const HubbardHamiltonian<double>& hamiltonian, const VectorX<double>& x) const;


    /*
     *  MARK: Dense unrestricted operator evaluations
     */

    /**
     *  Calculate the dense matrix representation of an unrestricted Hamiltonian in this ONV basis.
     *
     *  @param hamiltonian      An unrestricted Hamiltonian expressed in an orthonormal orbital basis.
     *
     *  @return A dense matrix represention of the Hamiltonian.
     */
    SquareMatrix<double> evaluateOperatorDense(const USQHamiltonian<double>& hamiltonian) const;


    /*
     *  MARK: Diagonal unrestricted operator evaluations
     */

    /**
     *  Calculate the diagonal of the dense matrix representation of an unrestricted Hamiltonian in this ONV basis.
     *
     *  @param hamiltonian      An unrestricted Hamiltonian expressed in an orthonormal orbital basis.
     *
     *  @return The diagonal of the dense matrix represention of the Hamiltonian.
     */
    VectorX<double> evaluateOperatorDiagonal(const USQHamiltonian<double>& hamiltonian) const;


    /*
     *  MARK: Unrestricted matrix-vector product evaluations
     */

    /**
     *  Calculate the matrix-vector product of (the matrix representation of) an unrestricted Hamiltonian with the given coefficient vector.
     *
     *  @param hamiltonian      An unrestricted Hamiltonian expressed in an orthonormal orbital basis.
     *  @param x                The coefficient vector of a linear expansion.
     *
     *  @return The coefficient vector of the linear expansion after being acted on with the given (matrix representation of) the Hamiltonian.
     */
    VectorX<double> evaluateOperatorMatrixVectorProduct(const USQHamiltonian<double>& usq_hamiltonian, const VectorX<double>& x) const;
};


}  // namespace GQCP
