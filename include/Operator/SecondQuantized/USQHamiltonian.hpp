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


#include "Basis/ScalarBasis.hpp"
#include "Basis/SingleParticleBasis.hpp"
#include "Basis/TransformationMatrix.hpp"
#include "HoppingMatrix.hpp"
#include "Molecule/Molecule.hpp"
#include "Operator/FirstQuantized/NuclearRepulsionOperator.hpp"
#include "Operator/FirstQuantized/OverlapOperator.hpp"
#include "Operator/SecondQuantized/SQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/SQTwoElectronOperator.hpp"
#include "RDM/OneRDM.hpp"
#include "RDM/TwoRDM.hpp"
#include "typedefs.hpp"
#include "Utilities/miscellaneous.hpp"


namespace GQCP {


/**
 *  A class for representing the second-quantized electronic Hamiltonian, it consists of one-electron and two-electron contributions
 *
 *  @tparam Scalar      the scalar type of the second-quantized parameters (i.e. the integrals)
 */
template <typename Scalar>
class USQHamiltonian {
private:
    SQHamiltonian<Scalar> sq_hamiltonian_alpha;  // alpha hamiltonian
    SQHamiltonian<Scalar> sq_hamiltonian_beta;  // beta hamiltonian

    ScalarSQTwoElectronOperator<Scalar> two_op_mixed;  // alpha & beta mixed operator represented as g_aabb

public:

    /*
     *  CONSTRUCTORS
     */
    USQHamiltonian() = default;

    /**
     *  @param sq_hamiltonian_alpha      alpha hamiltonian
     *  @param sq_hamiltonian_beta       beta hamiltonian
     *  @param two_op_mixed              alpha & beta mixed operator
     */
    USQHamiltonian(const SQHamiltonian<Scalar>& sq_hamiltonian_alpha, const SQHamiltonian<Scalar>& sq_hamiltonian_beta, const ScalarSQTwoElectronOperator<Scalar>& two_op_mixed) :
        sq_hamiltonian_alpha (sq_hamiltonian_alpha),
        sq_hamiltonian_beta (sq_hamiltonian_beta),
        two_op_mixed (two_op_mixed)
    {}

    /*
     *  NAMED CONSTRUCTORS
     */

    /**
     *  Construct the molecular Hamiltonian in a given single-particle basis
     *
     *  @param sp_basis     the single-particle basis in which the Hamiltonian should be expressed
     *
     *  @return a second-quantized molecular Hamiltonian. The molecular Hamiltonian has
     *      - one-electron contributions:
     *          - kinetic
     *          - nuclear attraction
     *      - two-electron contributions:
     *          - Coulomb repulsion
     *
     *  Note that this named constructor is only available for real matrix representations
     */
    template <typename Z = Scalar>
    static enable_if_t<std::is_same<Z, double>::value, USQHamiltonian<double>> Molecular(const SingleParticleBasis<Z, GTOShell>& sp_basis_alpha, const SingleParticleBasis<Z, GTOShell>& sp_basis_beta, const Molecule& molecule) {

        const SQHamiltonian<Scalar> sq_hamiltonian_alpha = SQHamiltonian<double>::Molecular(sp_basis_alpha, molecule);
        const SQHamiltonian<Scalar> sq_hamiltonian_beta = SQHamiltonian<double>::Molecular(sp_basis_beta, molecule);
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
     *  @return the contributions to the mixed alpha & beta two-electron part of the unrestricted Hamiltonian
     */
    const ScalarSQTwoElectronOperator<Scalar>& twoElectronMixed() const { return this->two_op_mixed; }

    /**
     *  In-place transform the matrix representations of Hamiltonian
     *
     *  @param T    the transformation matrix between the old and the new orbital basis
     */
    void transform(const TransformationMatrix<Scalar>& T) {

        this->sq_hamiltonian_alpha.transform(T);
        this->sq_hamiltonian_beta.transform(T);
        // Transform the mixed
        this->two_op_mixed.transform(T);
    }

    /**
     *  In-place transform the matrix representations of the alpha components of the unrestricted Hamiltonian
     *
     *  @param T    the transformation matrix between the old and the new orbital basis
     */
    void transformAlpha(const TransformationMatrix<Scalar>& T) {

        this->sq_hamiltonian_alpha.transform(T);
        // Transform the mixed
        auto new_two_electron_parameters = this->two_op_mixed.twoElectron().parameters();
        // transform the two electron parameters "g_aabb" to "g_a'a'bb"
        new_two_electron_parameters.template matrixContraction<Scalar>(T, 0);
        new_two_electron_parameters.template matrixContraction<Scalar>(T, 1);
        this->two_op_mixed = ScalarSQTwoElectronOperator<Scalar> ({new_two_electron_parameters});
    }

    /**
     *  In-place transform the matrix representations of the beta components of the unrestricted Hamiltonian
     *
     *  @param T    the transformation matrix between the old and the new orbital basis
     */
    void transformBeta(const TransformationMatrix<Scalar>& T) {

        this->sq_hamiltonian_beta.transform(T);
        // Transform the mixed
        auto new_two_electron_parameters = this->two_op_mixed.twoElectron().parameters();
         // transform the two electron parameters "g_aabb" to "g_aab'b'"
        new_two_electron_parameters.template matrixContraction<Scalar>(T, 2);
        new_two_electron_parameters.template matrixContraction<Scalar>(T, 3);
        this->two_op_mixed = ScalarSQTwoElectronOperator<Scalar> ({new_two_electron_parameters});
    }
};


}  // namespace GQCP
