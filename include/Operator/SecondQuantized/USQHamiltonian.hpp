// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2019  the GQCG developers
// 
// GQCG-gqcp is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// GQCG-gqcp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-gqcp.  If not, see <http://www.gnu.org/licenses/>.
// 
#pragma once


#include "Operator/SecondQuantized/SQHamiltonian.hpp"


namespace GQCP {


/**
 *  A class for representing the electronic Hamiltonian, expressed in an unrestricted spinor basis.
 * 
 *  @tparam Scalar      the scalar type of the second-quantized parameters (i.e. the integrals)
 */
template <typename Scalar>
class USQHamiltonian {
private:
    SQHamiltonian<Scalar> sq_hamiltonian_alpha;  // the alpha Hamiltonian
    SQHamiltonian<Scalar> sq_hamiltonian_beta;  // the beta Hamiltonian

    std::vector<ScalarSQTwoElectronOperator<Scalar>> two_op_mixed;  // the alpha & beta mixed two-electron operators (whose integrals are represented as g_aabb)
    ScalarSQTwoElectronOperator<Scalar> total_two_op_mixed;  // the total alpha & beta mixed two-electron operator (whose integrals are represented as g_aabb)

public:

    /*
     *  CONSTRUCTORS
     */
    USQHamiltonian() = default;

    /**
     *  @param sq_hamiltonian_alpha      the alpha Hamiltonian
     *  @param sq_hamiltonian_beta       the beta Hamiltonian
     *  @param two_op_mixed              the alpha & beta mixed two-electron operators (whose integrals are represented as g_aabb)
     */
    USQHamiltonian(const SQHamiltonian<Scalar>& sq_hamiltonian_alpha, const SQHamiltonian<Scalar>& sq_hamiltonian_beta, const std::vector<ScalarSQTwoElectronOperator<Scalar>>& two_op_mixed) :
        sq_hamiltonian_alpha (sq_hamiltonian_alpha),
        sq_hamiltonian_beta (sq_hamiltonian_beta),
        two_op_mixed (two_op_mixed)
    {
        // Check if the dimensions are compatible
        const auto dim = sq_hamiltonian_alpha.dimension();

        if (sq_hamiltonian_beta.dimension() != dim) {
            throw std::invalid_argument("USQHamiltonian::USQHamiltonian(const SQHamiltonian<Scalar>& sq_hamiltonian_alpha, const SQHamiltonian<Scalar>& sq_hamiltonian_beta, const ScalarSQTwoElectronOperator<Scalar>& two_op_mixed): The dimensions of the alpha and beta Hamiltonian are incompatible");
        }
        
        for (const auto& two_op : this->two_op_mixed) {
            if (two_op.dimension() != dim) {
                throw std::invalid_argument("USQHamiltonian::USQHamiltonian(const SQHamiltonian<Scalar>& sq_hamiltonian_alpha, const SQHamiltonian<Scalar>& sq_hamiltonian_beta, const ScalarSQTwoElectronOperator<Scalar>& two_op_mixed): The dimensions of the mixed two electron operator are incompatible with the Hamiltonian");
            }
        }
        
        // Calculate the total mixed two-electron operator
        this->total_two_op_mixed (dim);
        for (const auto& two_op : this->two_op_mixed) {
            total_two_op_mixed += two_op;
        }
    }

    /**
     *  @param sq_hamiltonian_alpha      the alpha Hamiltonian
     *  @param sq_hamiltonian_beta       the beta Hamiltonian
     *  @param two_op_mixed              the alpha & beta mixed two-electron operators (whose integrals are represented as g_aabb)
     */
    USQHamiltonian(const SQHamiltonian<Scalar>& sq_hamiltonian_alpha, const SQHamiltonian<Scalar>& sq_hamiltonian_beta, const ScalarSQTwoElectronOperator<Scalar>& two_op_mixed) :
        USQHamiltonian(sq_hamiltonian_alpha, sq_hamiltonian_beta, std::vector<ScalarSQTwoElectronOperator<Scalar>>({two_op_mixed}))
    {}


    /*
     *  NAMED CONSTRUCTORS
     */

    /**
     *  Construct the molecular Hamiltonian in a given single-particle basis
     *
     *  @param sp_basis           the initial single-particle basis in which both the alpha component and beta component of the Hamiltonian should be expressed
     *  @param molecule           the molecule on which the single particle is based
     *
     *  @return a second-quantized molecular unrestricted Hamiltonian. The molecular unrestricted Hamiltonian has:
     *      - restricted Hamiltonians for the alpha and beta components
     *      - mixed two-electron contributions
     *
     *  Note that this named constructor is only available for real matrix representations
     */
    template <typename Z = Scalar>
    static enable_if_t<std::is_same<Z, double>::value, USQHamiltonian<double>> Molecular(const SingleParticleBasis<Z, GTOShell>& sp_basis, const Molecule& molecule) {

        const SQHamiltonian<Scalar> sq_hamiltonian_alpha = SQHamiltonian<double>::Molecular(sp_basis, molecule);
        const SQHamiltonian<Scalar> sq_hamiltonian_beta = SQHamiltonian<double>::Molecular(sp_basis, molecule);

        // Initial basis for alpha and beta are identical so the mixed integrals are identical to spin specific components
        const ScalarSQTwoElectronOperator<double> two_op_mixed = sq_hamiltonian_alpha.twoElectron();

        return USQHamiltonian(sq_hamiltonian_alpha, sq_hamiltonian_beta, two_op_mixed);
    }

    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return the dimension of the Hamiltonian, i.e. the number of spinors in which it is expressed
     */
    size_t dimension() const { return this->sq_hamiltonian_alpha.dimension(); }

    /**
     *  @return the total of the alpha contributions of the unrestricted Hamiltonian
     */
    const SQHamiltonian<Scalar>& alphaHamiltonian() const { return this->sq_hamiltonian_alpha; }

    /**
     *  @return the total of the beta contributions of the unrestricted Hamiltonian
     */
    const SQHamiltonian<Scalar>& betaHamiltonian() const { return this->sq_hamiltonian_beta; }

    /**
     *  @return the total contributions to the mixed alpha & beta two-electron part of the unrestricted Hamiltonian
     */
    const ScalarSQTwoElectronOperator<Scalar>& twoElectronMixed() const { return this->total_two_op_mixed; }

    /**
     *  @return the total contributions to the mixed alpha & beta two-electron part of the unrestricted Hamiltonian
     */
    const std::vector<ScalarSQTwoElectronOperator<Scalar>>& twoElectronContributionsMixed() const { return this->two_op_mixed; }

    /**
     *  In-place transform the matrix representations of Hamiltonian
     *
     *  @param T    the transformation matrix between the old and the new orbital basis
     */
    void transform(const TransformationMatrix<Scalar>& T) {

        this->sq_hamiltonian_alpha.transform(T);
        this->sq_hamiltonian_beta.transform(T);

        // Transform the mixed
        this->total_two_op_mixed.transform(T);
        for (auto& two_op : this->two_op_mixed) {
            two_op.transform(T);
        }
    }

    /**
     *  In-place transform the matrix representations of the alpha components of the unrestricted Hamiltonian
     *
     *  @param T    the transformation matrix between the old and the new orbital basis
     */
    void transformAlpha(const TransformationMatrix<Scalar>& T) {

        this->sq_hamiltonian_alpha.transform(T);

        // Transform the mixed two-electron operators and their total
        auto new_two_electron_parameters = this->total_two_op_mixed.parameters();

        // transform the two electron parameters "g_aabb" to "g_a'a'bb"
        new_two_electron_parameters.template matrixContraction<Scalar>(T, 0);
        new_two_electron_parameters.template matrixContraction<Scalar>(T, 1);
        this->total_two_op_mixed = ScalarSQTwoElectronOperator<Scalar> ({new_two_electron_parameters});

        for (auto& two_op : this->two_op_mixed) {

            auto new_two_electron_parameters = two_op.parameters();
            // transform the two electron parameters "g_aabb" to "g_a'a'bb"
            new_two_electron_parameters.template matrixContraction<Scalar>(T, 0);
            new_two_electron_parameters.template matrixContraction<Scalar>(T, 1);
            two_op = ScalarSQTwoElectronOperator<Scalar> ({new_two_electron_parameters});
        }
    }


    /**
     *  In-place transform the matrix representations of the beta components of the unrestricted Hamiltonian
     *
     *  @param T    the transformation matrix between the old and the new orbital basis
     */
    void transformBeta(const TransformationMatrix<Scalar>& T) {

        this->sq_hamiltonian_beta.transform(T);

        // Transform the mixed two-electron operators and their total
        auto new_two_electron_parameters = this->total_two_op_mixed.parameters();

        // transform the two electron parameters "g_aabb" to "g_aab'b'"
        new_two_electron_parameters.template matrixContraction<Scalar>(T, 2);
        new_two_electron_parameters.template matrixContraction<Scalar>(T, 3);
        this->total_two_op_mixed = ScalarSQTwoElectronOperator<Scalar> ({new_two_electron_parameters});

        for (auto& two_op : this->two_op_mixed) {

            auto new_two_electron_parameters = two_op.parameters();
            // transform the two electron parameters "g_aabb" to "g_aab'b'"
            new_two_electron_parameters.template matrixContraction<Scalar>(T, 2);
            new_two_electron_parameters.template matrixContraction<Scalar>(T, 3);
            two_op = ScalarSQTwoElectronOperator<Scalar> ({new_two_electron_parameters});
        }
    }


    /**
     *  Constrain the alpha component Hamiltonian according to the convention: - lambda * constraint
     *
     *  @param one_op   the one-electron operator used as a constraint
     *  @param lambda   the Lagrangian multiplier for the constraint
     *
     *  @return a copy of the constrained Hamiltonian
     *
     *  Note that this method is only available for real matrix representations
     */
    template<typename Z = Scalar>
    enable_if_t<std::is_same<Z, double>::value, USQHamiltonian<double>> constrainAlpha(const ScalarSQOneElectronOperator<double>& one_op, double lambda) const {

        auto const constrained_component = this->sq_hamiltonian_alpha.constrain(one_op, lambda);

        return USQHamiltonian(constrained_component, this->sq_hamiltonian_beta, this->two_op_mixed);
    }

    /**
     *  Constrain the beta component Hamiltonian according to the convention: - lambda * constraint
     *
     *  @param one_op   the one-electron operator used as a constraint
     *  @param lambda   The Lagrangian multiplier for the constraint
     *
     *  @return a copy of the constrained Hamiltonian
     *
     *  Note that this method is only available for real matrix representations
     */
    template<typename Z = Scalar>
    enable_if_t<std::is_same<Z, double>::value, USQHamiltonian<double>> constrainBeta(const ScalarSQOneElectronOperator<double>& one_op, double lambda) const {

        auto const constrained_component = this->sq_hamiltonian_beta.constrain(one_op, lambda);

        return USQHamiltonian(this->sq_hamiltonian_alpha, constrained_component, this->two_op_mixed);
    }
};


}  // namespace GQCP
