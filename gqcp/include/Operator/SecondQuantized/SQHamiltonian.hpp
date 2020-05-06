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


#include "Basis/ScalarBasis/ScalarBasis.hpp"
#include "Basis/SpinorBasis/JacobiRotationParameters.hpp"
#include "Basis/SpinorBasis/RSpinorBasis.hpp"
#include "Basis/TransformationMatrix.hpp"
#include "Molecule/Molecule.hpp"
#include "Operator/FirstQuantized/NuclearRepulsionOperator.hpp"
#include "Operator/FirstQuantized/OverlapOperator.hpp"
#include "Operator/SecondQuantized/ModelHamiltonian/HubbardHamiltonian.hpp"
#include "Operator/SecondQuantized/SQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/SQTwoElectronOperator.hpp"
#include "Processing/RDM/OneRDM.hpp"
#include "Processing/RDM/TwoRDM.hpp"
#include "Utilities/miscellaneous.hpp"
#include "Utilities/type_traits.hpp"


namespace GQCP {


/**
 *  A class for representing the second-quantized electronic Hamiltonian, it consists of one-electron and two-electron contributions
 *
 *  @tparam Scalar      the scalar type of the second-quantized parameters (i.e. the integrals)
 */
template <typename Scalar>
class SQHamiltonian {
private:
    ScalarSQOneElectronOperator<Scalar> total_one_op;  // one-electron interactions (i.e. the core Hamiltonian)
    ScalarSQTwoElectronOperator<Scalar> total_two_op;  // two-electron interactions

    std::vector<ScalarSQOneElectronOperator<Scalar>> one_ops;  // the core (i.e. one-electron) contributions to the Hamiltonian
    std::vector<ScalarSQTwoElectronOperator<Scalar>> two_ops;  // the two-electron contributions to the Hamiltonian


public:
    /*
     *  CONSTRUCTORS
     */
    SQHamiltonian() = default;

    /**
     *  @param one_ops      the core (i.e. one-electron) contributions to the Hamiltonian
     *  @param two_ops      the two-electron contributions to the Hamiltonian
     */
    SQHamiltonian(const std::vector<ScalarSQOneElectronOperator<Scalar>>& one_ops, const std::vector<ScalarSQTwoElectronOperator<Scalar>>& two_ops) :
        one_ops {one_ops},
        two_ops {two_ops} {

        // Check if the dimensions are compatible
        const std::invalid_argument dimension_error("SQHamiltonian::SQHamiltonian(const std::vector<ScalarSQOneElectronOperator<Scalar>& one_ops, const std::vector<ScalarSQTwoElectronOperator<Scalar>& two_ops: The dimensions of the operators and coefficients matrix are incompatible");

        const auto dim = one_ops[0].dimension();
        for (const auto& one_op : this->one_ops) {
            if (one_op.dimension() != dim) {
                throw dimension_error;
            }
        }

        for (const auto& two_op : this->two_ops) {
            if (two_op.dimension() != dim) {
                throw dimension_error;
            }
        }


        // Calculate the total one-electron operator
        QCMatrix<Scalar> total_one_op_par {dim};
        total_one_op_par.setZero();
        for (const auto& one_op : this->one_ops) {
            total_one_op_par += one_op.parameters();
        }
        this->total_one_op = ScalarSQOneElectronOperator<Scalar>(total_one_op_par);


        // Calculate the total two-electron operator
        QCRankFourTensor<Scalar> total_two_op_par(dim);
        total_two_op_par.setZero();
        for (const auto& two_op : this->two_ops) {
            total_two_op_par += two_op.parameters().Eigen();
        }
        this->total_two_op = ScalarSQTwoElectronOperator<Scalar>(total_two_op_par);
    }


    /**
     *  @param h            the (total) one-electron (i.e. core) integrals
     *  @param g            the (total) two-electron integrals
     */
    SQHamiltonian(const ScalarSQOneElectronOperator<Scalar>& h, const ScalarSQTwoElectronOperator<Scalar>& g) :
        SQHamiltonian(std::vector<ScalarSQOneElectronOperator<Scalar>>{h}, std::vector<ScalarSQTwoElectronOperator<Scalar>>{g}) {}


    /*
     *  NAMED CONSTRUCTORS
     */

    /**
     *  @param hubbard_hamiltonian              a Hubbard model Hamiltonian
     *
     *  @return a full SQHamiltonian generated from the given Hubbard model Hamiltonian
     *
     *  @note This named constructor is only available for real matrix representations.
     */
    template <typename Z = Scalar>
    static enable_if_t<std::is_same<Z, double>::value, SQHamiltonian<double>> FromHubbard(const GQCP::HubbardHamiltonian<double>& hubbard_hamiltonian) {

        const auto h = hubbard_hamiltonian.core();
        const auto g = hubbard_hamiltonian.twoElectron();

        return SQHamiltonian(h, g);
    }


    /**
     *  Construct the molecular Hamiltonian in a given spinor basis.
     *
     *  @param spinor_basis     the spinor basis in which the Hamiltonian should be expressed
     *  @param molecule         the molecule on which the single particle is based
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
    static enable_if_t<std::is_same<Z, double>::value, SQHamiltonian<double>> Molecular(const RSpinorBasis<Z, GTOShell>& spinor_basis, const Molecule& molecule) {

        // Calculate the integrals for the molecular Hamiltonian
        const auto T = spinor_basis.quantize(Operator::Kinetic());
        const auto V = spinor_basis.quantize(Operator::NuclearAttraction(molecule));
        ScalarSQOneElectronOperator<double> H = T + V;

        const auto g = spinor_basis.quantize(Operator::Coulomb());

        return SQHamiltonian(H, g);
    }


    /**
     *  @param K        the number of orbitals
     *
     *  @return a random Hamiltonian with values uniformly distributed between [-1,1]
     *
     *  Note that this named constructor is only available for real representations
     */
    template <typename Z = Scalar>
    static enable_if_t<std::is_same<Z, double>::value, SQHamiltonian<double>> Random(size_t K) {

        ScalarSQOneElectronOperator<double> H {QCMatrix<double>::Random(K, K)};  // uniformly distributed between [-1,1]


        // Unfortunately, the Tensor module provides uniform random distributions between [0, 1]
        QCRankFourTensor<double> g {K};
        g.setRandom();

        // Move the distribution from [0, 1] -> [-1, 1]
        for (size_t i = 0; i < K; i++) {
            for (size_t j = 0; j < K; j++) {
                for (size_t k = 0; k < K; k++) {
                    for (size_t l = 0; l < K; l++) {
                        g(i, j, k, l) = 2 * g(i, j, k, l) - 1;  // scale from [0, 1] -> [0, 2] -> [-1, 1]
                    }
                }
            }
        }

        return SQHamiltonian<double>(H, ScalarSQTwoElectronOperator<double>(g));
    }


    /**
     *  @param fcidump_file     the name of the FCIDUMP file
     *
     *  @return the Hamiltonian corresponding to the contents of an FCIDUMP file
     *
     *  Note that this named constructor is only available for real matrix representations
     */
    template <typename Z = Scalar>
    static enable_if_t<std::is_same<Z, double>::value, SQHamiltonian<double>> ReadFCIDUMP(const std::string& fcidump_file) {

        std::ifstream input_file_stream = validateAndOpen(fcidump_file, "FCIDUMP");


        // Do the actual parsing

        //  Get the number of orbitals to check if it's a valid FCIDUMP file
        std::string start_line;  // first line contains orbitals and electron count
        std::getline(input_file_stream, start_line);
        std::stringstream linestream {start_line};

        size_t K = 0;
        char iter;

        while (linestream >> iter) {
            if (iter == '=') {
                linestream >> K;  // right here we have the number of orbitals
                break;            // we can finish reading the linestream after we found K
            }
        }

        if (K == 0) {
            throw std::invalid_argument("SQHamiltonian::ReadFCIDUMP(std::string): The .FCIDUMP-file is invalid: could not read a number of orbitals.");
        }


        QCMatrix<double> h_core = QCMatrix<double>::Zero(K, K);
        QCRankFourTensor<double> g {K};
        g.setZero();

        //  Skip 3 lines
        for (size_t counter = 0; counter < 3; counter++) {
            std::getline(input_file_stream, start_line);
        }


        //  Start reading in the one- and two-electron integrals
        double x;
        size_t i, j, a, b;

        std::string line;
        while (std::getline(input_file_stream, line)) {
            std::istringstream iss {line};

            // Based on what the values of the indices are, we can read one-electron integrals, two-electron integrals and the internuclear repulsion energy
            //  See also (http://hande.readthedocs.io/en/latest/manual/integrals.html)
            //  I think the documentation is a bit unclear for the two-electron integrals, but we can rest assured that FCIDUMP files give the two-electron integrals in CHEMIST's notation.
            iss >> x >> i >> a >> j >> b;

            //  Single-particle eigenvalues (skipped)
            if ((a == 0) && (j == 0) && (b == 0)) {
            }

            //  One-electron integrals (h_core)
            else if ((j == 0) && (b == 0)) {
                size_t p = i - 1;
                size_t q = a - 1;
                h_core(p, q) = x;

                // Apply the permutational symmetry for real orbitals
                h_core(q, p) = x;
            }

            //  Two-electron integrals are given in CHEMIST'S NOTATION, so just copy them over
            else if ((i > 0) && (a > 0) && (j > 0) && (b > 0)) {
                size_t p = i - 1;
                size_t q = a - 1;
                size_t r = j - 1;
                size_t s = b - 1;
                g(p, q, r, s) = x;

                // Apply the permutational symmetries for real orbitals
                g(p, q, s, r) = x;
                g(q, p, r, s) = x;
                g(q, p, s, r) = x;

                g(r, s, p, q) = x;
                g(s, r, p, q) = x;
                g(r, s, q, p) = x;
                g(s, r, q, p) = x;
            }
        }  // while loop


        return SQHamiltonian(ScalarSQOneElectronOperator<Scalar>(h_core), ScalarSQTwoElectronOperator<Scalar>(g));
    }


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @param N_P      the number of electron pairs
     *
     *  @return the Edmiston-Ruedenberg localization index g(i,i,i,i)
     *
     *  Note that this method is only available for real matrix representations
     */
    template <typename Z = Scalar>
    enable_if_t<std::is_same<Z, double>::value, double> calculateEdmistonRuedenbergLocalizationIndex(size_t N_P) const {

        const auto& g_par = this->total_two_op.parameters();

        // TODO: when Eigen releases TensorTrace, use it here
        double localization_index = 0.0;
        for (size_t i = 0; i < N_P; i++) {
            localization_index += g_par(i, i, i, i);
        }

        return localization_index;
    }


    /**
     *  @return the effective one-electron integrals
     */
    ScalarSQOneElectronOperator<Scalar> calculateEffectiveOneElectronIntegrals() const {

        return this->core() + this->twoElectron().effectiveOneElectronPartition();
    }


    /**
     *  @param D            the 1-RDM
     *  @param d            the 2-RDM
     *
     *  @return the expectation value of this Hamiltonian
     */
    Scalar calculateExpectationValue(const OneRDM<Scalar>& D, const TwoRDM<Scalar>& d) const {

        return this->core().calculateExpectationValue(D)[0] + this->twoElectron().calculateExpectationValue(d)[0];  // SQHamiltonian contains ScalarSQOperators, so we access with [0]
    }


    /**
     *  @param D      the 1-DM (or the response 1-DM for made-variational wave function models)
     *  @param d      the 2-DM (or the response 2-DM for made-variational wave function models)
     *
     *  @return the (generalized) Fockian matrix
     */
    SquareMatrix<Scalar> calculateFockianMatrix(const OneRDM<double>& D, const TwoRDM<double>& d) const {

        return this->core().calculateFockianMatrix(D, d)[0] + this->twoElectron().calculateFockianMatrix(D, d)[0];  // SQHamiltonian has one- and two-electron contributions, so access with [0] accordingly
    }


    /**
     *  @param N_P          the number of electron pairs
     * 
     *  @return the inactive Fockian matrix
     */
    ScalarSQOneElectronOperator<Scalar> calculateInactiveFockian(const size_t N_P) const {

        const auto& h_par = this->core().parameters();
        const auto& g_par = this->twoElectron().parameters();


        // A KISS implementation of the calculation of the inactive Fockian matrix
        auto F_par = h_par;  // one-electron part

        // Two-electron part
        for (size_t p = 0; p < this->dimension(); p++) {
            for (size_t q = 0; q < this->dimension(); q++) {

                for (size_t i = 0; i < N_P; i++) {
                    F_par(p, q) += 2 * g_par(p, q, i, i) - g_par(p, i, i, q);
                }
            }
        }  // F elements loop

        return ScalarSQOneElectronOperator<Scalar>(F_par);
    }


    /**
     *  @param D      the 1-DM (or the response 1-DM for made-variational wave function models)
     *  @param d      the 2-DM (or the response 2-DM for made-variational wave function models)
     *
     *  @return the (generalized) super-Fockian matrix
     */
    SquareRankFourTensor<Scalar> calculateSuperFockianMatrix(const OneRDM<double>& D, const TwoRDM<double>& d) const {

        return this->core().calculateSuperFockianMatrix(D, d)[0].Eigen() + this->twoElectron().calculateSuperFockianMatrix(D, d)[0].Eigen();  // SQHamiltonian contains ScalarSQOperators
    }


    /**
     *  @return the 'core' Hamiltonian, i.e. the total of the one-electron contributions to the Hamiltonian
     */
    const ScalarSQOneElectronOperator<Scalar>& core() const { return this->total_one_op; }

    /**
     *  @return the contributions to the 'core' Hamiltonian
     */
    const std::vector<ScalarSQOneElectronOperator<Scalar>>& coreContributions() const { return this->one_ops; }

    /**
     *  @return the dimension of the Hamiltonian, i.e. the number of spinors in which it is expressed
     */
    size_t dimension() const { return this->core().dimension(); }

    /**
     *  Using a random rotation matrix, transform the matrix representations of the Hamiltonian
     *
     *  Note that this method is only available for real matrix representations
     */
    template <typename Z = Scalar>
    enable_if_t<std::is_same<Z, double>::value> randomRotate() {

        // Get a random unitary matrix by diagonalizing a random symmetric matrix
        const auto K = this->dimension();
        TransformationMatrix<double> A_random = TransformationMatrix<double>::Random(K, K);
        TransformationMatrix<double> A_symmetric = A_random + A_random.transpose();
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> unitary_solver {A_symmetric};
        TransformationMatrix<double> U_random = unitary_solver.eigenvectors();

        this->rotate(U_random);
    }


    /**
     *  In-place rotate the matrix representations of Hamiltonian
     *
     *  @param U    the unitary rotation matrix between the old and the new orbital basis
     */
    void rotate(const TransformationMatrix<Scalar>& U) {

        // Rotate the one-electron contributions
        for (auto& one_op : this->one_ops) {
            one_op.rotate(U);
        }

        // Rotate the two-electron contributions
        for (auto& two_op : this->two_ops) {
            two_op.rotate(U);
        }

        // Rotate the totals
        this->total_one_op.rotate(U);
        this->total_two_op.rotate(U);
    }


    /**
     *  In-place rotate the matrix representations of the Hamiltonian using a unitary Jacobi rotation matrix constructed from the Jacobi rotation parameters. Note that this function is only available for real (double) matrix representations
     *
     *  @param jacobi_rotation_parameters       the Jacobi rotation parameters (p, q, angle) that are used to specify a Jacobi rotation: we use the (cos, sin, -sin, cos) definition for the Jacobi rotation matrix
     */
    template <typename Z = Scalar>
    enable_if_t<std::is_same<Z, double>::value> rotate(const JacobiRotationParameters& jacobi_rotation_parameters) {

        // Transform the one-electron contributions
        for (auto& one_op : this->one_ops) {
            one_op.rotate(jacobi_rotation_parameters);
        }

        // Transform the two-electron contributions
        for (auto& two_op : this->two_ops) {
            two_op.rotate(jacobi_rotation_parameters);
        }

        // Transform the totals
        this->total_one_op.rotate(jacobi_rotation_parameters);
        this->total_two_op.rotate(jacobi_rotation_parameters);
    }


    /**
     *  In-place transform the matrix representations of Hamiltonian
     *
     *  @param T    the transformation matrix between the old and the new orbital basis
     */
    void transform(const TransformationMatrix<Scalar>& T) {

        // Transform the one-electron contributions
        for (auto& one_op : this->one_ops) {
            one_op.transform(T);
        }

        // Transform the two-electron contributions
        for (auto& two_op : this->two_ops) {
            two_op.transform(T);
        }

        // Transform the totals
        this->total_one_op.transform(T);
        this->total_two_op.transform(T);
    }


    /**
     *  @return the total of the two-electron contributions to the Hamiltonian
     */
    const ScalarSQTwoElectronOperator<Scalar>& twoElectron() const { return this->total_two_op; }

    /**
     *  @return the contributions to the two-electron part of the Hamiltonian
     */
    const std::vector<ScalarSQTwoElectronOperator<Scalar>>& twoElectronContributions() const { return this->two_ops; }
};


/*
 *  OPERATORS
 */

/**
 *  Add the second-quantized forms of a (scalar) one-electron operator to that of a Hamiltonian.
 * 
 *  @tparam Scalar              the type that is used to represent elements of the (scalar) one-electron operator and the Hamiltonian
 * 
 *  @param sq_hamiltonian       the second-quantized Hamiltonian
 *  @param sq_one_op            the (scalar) second-quantized one-electron operator
 * 
 *  @return a new second-quantized Hamiltonian
 */
template <typename Scalar>
SQHamiltonian<Scalar> operator+(const SQHamiltonian<Scalar>& sq_hamiltonian, const ScalarSQOneElectronOperator<Scalar>& sq_one_op) {

    // Make a copy of the one-electron part in order to create a new Hamiltonian
    auto sq_one_ops = sq_hamiltonian.coreContributions();

    // 'Add' the one-electron operator
    sq_one_ops.push_back(sq_one_op);

    return SQHamiltonian<Scalar>(sq_one_ops, sq_hamiltonian.twoElectronContributions());
}


/**
 *  Subtract a (scalar) second-quantized one-electron operator from a second-quantized Hamiltonian.
 * 
 *  @tparam Scalar              the type that is used to represent elements of the (scalar) one-electron operator and the Hamiltonian
 * 
 *  @param sq_hamiltonian       the second-quantized Hamiltonian
 *  @param sq_one_op            the (scalar) second-quantized one-electron operator
 * 
 *  @return a new second-quantized Hamiltonian
 */
template <typename Scalar>
SQHamiltonian<Scalar> operator-(const SQHamiltonian<Scalar>& sq_hamiltonian, const ScalarSQOneElectronOperator<Scalar>& sq_one_op) {

    return sq_hamiltonian + (-sq_one_op);
}


/**
 *  Add the second-quantized forms of a (scalar) two-electron operator to that of a Hamiltonian.
 * 
 *  @tparam Scalar              the type that is used to represent elements of the (scalar) two-electron operator and the Hamiltonian
 * 
 *  @param sq_hamiltonian       the second-quantized Hamiltonian
 *  @param sq_two_op            the (scalar) second-quantized two-electron operator
 * 
 *  @return a new second-quantized Hamiltonian
 */
template <typename Scalar>
SQHamiltonian<Scalar> operator+(const SQHamiltonian<Scalar>& sq_hamiltonian, const ScalarSQTwoElectronOperator<Scalar>& sq_two_op) {

    // Make a copy of the two-electron part in order to create a new Hamiltonian
    auto sq_two_ops = sq_hamiltonian.twoElectronContributions();

    // 'Add' the two-electron operator
    sq_two_ops.push_back(sq_two_ops);

    return SQHamiltonian<Scalar>(sq_hamiltonian.coreContributions(), sq_two_ops);
}


/**
 *  Subtract a (scalar) second-quantized two-electron operator from a second-quantized Hamiltonian.
 * 
 *  @tparam Scalar              the type that is used to represent elements of the (scalar) two-electron operator and the Hamiltonian
 * 
 *  @param sq_hamiltonian       the second-quantized Hamiltonian
 *  @param sq_two_op            the (scalar) second-quantized two-electron operator
 * 
 *  @return a new second-quantized Hamiltonian
 */
template <typename Scalar>
SQHamiltonian<Scalar> operator-(const SQHamiltonian<Scalar>& sq_hamiltonian, const ScalarSQTwoElectronOperator<Scalar>& sq_two_op) {

    return sq_hamiltonian + (-sq_two_op);
}


}  // namespace GQCP
