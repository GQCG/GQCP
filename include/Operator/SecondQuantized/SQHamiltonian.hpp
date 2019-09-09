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
#include "OrbitalOptimization/JacobiRotationParameters.hpp"
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
class SQHamiltonian {
private:
    size_t K;  // the number of spatial orbitals

    std::shared_ptr<ScalarBasis<GTOShell>> ao_basis;  // the initial atomic orbitals
    TransformationMatrix<Scalar> T_total;  // total transformation matrix between the current (restricted) molecular orbitals and the atomic orbitals
    ScalarSQOneElectronOperator<Scalar> S;  // overlap

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
     *  @param ao_basis     the initial AO basis
     *  @param S            the overlap integrals
     *  @param one_ops      the core (i.e. one-electron) contributions to the Hamiltonian
     *  @param g            the two-electron contributions to the Hamiltonian
     *  @param T            the transformation matrix between the current molecular orbitals and the atomic orbitals
     */
    SQHamiltonian(std::shared_ptr<ScalarBasis<GTOShell>> ao_basis, const ScalarSQOneElectronOperator<Scalar>& S, const std::vector<ScalarSQOneElectronOperator<Scalar>>& one_ops, const std::vector<ScalarSQTwoElectronOperator<Scalar>>& two_ops, const TransformationMatrix<Scalar>& T) :
        ao_basis (std::move(ao_basis)),
        K (S.dimension()),
        S (S),
        one_ops (one_ops),
        two_ops (two_ops),
        T_total (T)
    {
        const std::string function_intro ("SQHamiltonian::SQHamiltonian(std::shared_ptr<ScalarBasis<GTOShell>> ao_basis, const ScalarSQOneElectronOperator<Scalar>& S, const std::vector<ScalarSQOneElectronOperator<Scalar>& one_ops, const std::vector<ScalarSQTwoElectronOperator<Scalar>& two_ops, const TransformationMatrix<Scalar>& T): ");

        // Check if the dimensions are compatible
        const std::invalid_argument dimension_error (function_intro + "The dimensions of the operators and coefficient matrix are incompatible.");

        if (this->ao_basis) {  // ao_basis is not nullptr
            if (this->dimension() != this->ao_basis->numberOfBasisFunctions()) {
                throw dimension_error;
            }
        }

        const auto dimension_of_first = one_ops[0].dimension();
        for (const auto& one_op : this->one_ops) {
            if (one_op.dimension() != dimension_of_first) {
                throw dimension_error;
            }
        }

        for (const auto& two_op : this->two_ops) {
            if (two_op.dimension() != dimension_of_first) {
                throw dimension_error;
            }
        }

        if ((T.cols() != this->dimension()) || (T.rows() != this->dimension())) {
            throw dimension_error;
        }


        // Calculate the total one-electron operator
        QCMatrix<Scalar> total_one_op_par (this->dimension());
        total_one_op_par.setZero();
        for (const auto& one_op : this->one_ops) {
            total_one_op_par += one_op.parameters();
        }
        this->total_one_op = ScalarSQOneElectronOperator<Scalar>({total_one_op_par});


        // Calculate the total two-electron operator
        QCRankFourTensor<Scalar> total_two_op_par (this->dimension());
        total_two_op_par.setZero();
        for (const auto& two_op : this->two_ops) {
            total_two_op_par += two_op.parameters().Eigen();
        }
        this->total_two_op = ScalarSQTwoElectronOperator<Scalar>({total_two_op_par});



        // Check if the underlying overlap matrix is not a zero matrix
        if (S.parameters().isZero(1.0e-08)) {
            throw std::invalid_argument(function_intro + "The underlying overlap matrix cannot be a zero matrix.");
        }
    }


    /**
     *  @param ao_basis     the initial AO basis
     *  @param S            the overlap integrals
     *  @param h            the one-electron integrals H_core
     *  @param g            the two-electron integrals
     *  @param T            the transformation matrix between the current molecular orbitals and the atomic orbitals
     */
    SQHamiltonian(std::shared_ptr<ScalarBasis<GTOShell>> ao_basis, const ScalarSQOneElectronOperator<Scalar>& S, const ScalarSQOneElectronOperator<Scalar>& h, const ScalarSQTwoElectronOperator<Scalar>& g, const TransformationMatrix<Scalar>& T) :
        SQHamiltonian(ao_basis, S, std::vector<ScalarSQOneElectronOperator<Scalar>>({h}), std::vector<ScalarSQTwoElectronOperator<Scalar>>({g}), T)
    {}


    /**
     *  A constructor that transforms the given Hamiltonian with a transformation matrix
     *
     *  @param sq_hamiltonian       the current Hamiltonian
     *  @param C                    the transformation matrix to be applied to the given Hamiltonian
     */
    SQHamiltonian(const SQHamiltonian<Scalar>& sq_hamiltonian, const TransformationMatrix<Scalar>& C) :
        SQHamiltonian<Scalar>(sq_hamiltonian.ao_basis, sq_hamiltonian.S, sq_hamiltonian.core(), sq_hamiltonian.twoElectron(), sq_hamiltonian.T_total)
    {
        // We have now initialized the new Hamiltonian to be a copy of the given Hamiltonian, so now we will transform
        this->transform(C);
    }



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
    static enable_if_t<std::is_same<Z, double>::value, SQHamiltonian<double>> Molecular(const SingleParticleBasis<Z, GTOShell>& sp_basis, const Molecule& molecule) {

        // Calculate the integrals for the molecular Hamiltonian
        const auto S = sp_basis.quantize(Operator::Overlap());
        const auto T = sp_basis.quantize(Operator::Kinetic());
        const auto V = sp_basis.quantize(Operator::NuclearAttraction(molecule));
        ScalarSQOneElectronOperator<double> H = T + V;

        const auto g = sp_basis.quantize(Operator::Coulomb());


        // Construct the initial transformation matrix: the identity matrix
        auto nbf = sp_basis.numberOfBasisFunctions();
        TransformationMatrix<double> T_total = TransformationMatrix<double>::Identity(nbf, nbf);

        return SQHamiltonian(std::make_shared<ScalarBasis<GTOShell>>(sp_basis.scalarBasis()), S, H, g, T_total);
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

        ScalarSQOneElectronOperator<double> S ({QCMatrix<double>::Identity(K, K)});  // the underlying orbital basis can be chosen as orthonormal, since the form of the underlying orbitals doesn't really matter
        TransformationMatrix<double> C = TransformationMatrix<double>::Identity(K, K);  // the transformation matrix here doesn't really mean anything, because it doesn't link to any AO basis

        ScalarSQOneElectronOperator<double> H ({QCMatrix<double>::Random(K, K)});  // uniformly distributed between [-1,1]


        // Unfortunately, the Tensor module provides uniform random distributions between [0, 1]
        QCRankFourTensor<double> g (K);
        g.setRandom();

        // Move the distribution from [0, 1] -> [-1, 1]
        for (size_t i = 0; i < K; i++) {
            for (size_t j = 0; j < K; j++) {
                for (size_t k = 0; k < K; k++) {
                    for (size_t l = 0; l < K; l++) {
                        g(i,j,k,l) = 2*g(i,j,k,l) - 1;  // scale from [0, 1] -> [0, 2] -> [-1, 1]
                    }
                }
            }
        }

        std::shared_ptr<ScalarBasis<GTOShell>> ao_basis;  // nullptr because it doesn't make sense to set a scalar basis

        return SQHamiltonian<double>(ao_basis, S, H, ScalarSQTwoElectronOperator<double>({g}), C);
    }


    /**
     *  @param fcidump_file     the name of the FCIDUMP file
     *
     *  @return the Hamiltonian corresponding to the contents of an FCIDUMP file
     *
     *  Note that this named constructor is only available for real matrix representations
     */
    template<typename Z = Scalar>
    static enable_if_t<std::is_same<Z, double>::value, SQHamiltonian<double>> ReadFCIDUMP(const std::string& fcidump_file) {

        std::ifstream input_file_stream = validateAndOpen(fcidump_file, "FCIDUMP");


        // Do the actual parsing

        //  Get the number of orbitals to check if it's a valid FCIDUMP file
        std::string start_line;  // first line contains orbitals and electron count
        std::getline(input_file_stream, start_line);
        std::stringstream linestream (start_line);

        size_t K = 0;
        char iter;

        while (linestream >> iter) {
            if (iter == '=') {
                linestream >> K;  // right here we have the number of orbitals
                break;  // we can finish reading the linestream after we found K
            }
        }

        if (K == 0) {
            throw std::invalid_argument("SQHamiltonian::ReadFCIDUMP(std::string): The .FCIDUMP-file is invalid: could not read a number of orbitals.");
        }


        QCMatrix<double> h_core = QCMatrix<double>::Zero(K, K);
        QCRankFourTensor<double> g (K);
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
            std::istringstream iss (line);

            // Based on what the values of the indices are, we can read one-electron integrals, two-electron integrals and the internuclear repulsion energy
            //  See also (http://hande.readthedocs.io/en/latest/manual/integrals.html)
            //  I think the documentation is a bit unclear for the two-electron integrals, but we can rest assured that FCIDUMP files give the two-electron integrals in CHEMIST's notation.
            iss >> x >> i >> a >> j >> b;

            //  Single-particle eigenvalues (skipped)
            if ((a == 0) && (j == 0) && (b == 0)) {}

            //  One-electron integrals (h_core)
            else if ((j == 0) && (b == 0)) {
                size_t p = i - 1;
                size_t q = a - 1;
                h_core(p,q) = x;

                // Apply the permutational symmetry for real orbitals
                h_core(q,p) = x;
            }

            //  Two-electron integrals are given in CHEMIST'S NOTATION, so just copy them over
            else if ((i > 0) && (a > 0) && (j > 0) && (b > 0)) {
                size_t p = i - 1;
                size_t q = a - 1;
                size_t r = j - 1;
                size_t s = b - 1;
                g(p,q,r,s) = x;

                // Apply the permutational symmetries for real orbitals
                g(p,q,s,r) = x;
                g(q,p,r,s) = x;
                g(q,p,s,r) = x;

                g(r,s,p,q) = x;
                g(s,r,p,q) = x;
                g(r,s,q,p) = x;
                g(s,r,q,p) = x;
            }
        }  // while loop


        // Make the ingredients to construct SQHamiltonian
        std::shared_ptr<ScalarBasis<GTOShell>> ao_basis;  // nullptr
        ScalarSQOneElectronOperator<Scalar> S ({QCMatrix<double>::Identity(K, K)});
        TransformationMatrix<double> C = TransformationMatrix<double>::Identity(K, K);

        return SQHamiltonian(ao_basis, S, ScalarSQOneElectronOperator<Scalar>({h_core}), ScalarSQTwoElectronOperator<Scalar>({g}), C);
    }


    /**
     *  @param H      a Hubbard hopping matrix
     *
     *  @return a Hubbard Hamiltonian generated from the Hubbard hopping matrix
     *
     *  Note that this named constructor is only available for real matrix representations
     */
    template<typename Z = Scalar>
    static enable_if_t<std::is_same<Z, double>::value, SQHamiltonian<double>> Hubbard(const HoppingMatrix& H) {

        const size_t K = H.numberOfLatticeSites();

        QCMatrix<double> h = QCMatrix<double>::Zero(K, K);
        QCRankFourTensor<double> g (K);
        g.setZero();


        for (size_t i = 0; i < K; i++) {
            for (size_t j = i; j < K; j++) {
                if (i == j) {
                    g(i,i,i,i) = H(i,i);
                } else {
                    h(i,j) = H(i,j);
                    h(j,i) = H(j,i);
                }
            }
        }


        // Make the ingredients to construct SQHamiltonian
        std::shared_ptr<ScalarBasis<GTOShell>> ao_basis;  // nullptr
        ScalarSQOneElectronOperator<double> S ({QCMatrix<double>::Identity(K, K)});
        TransformationMatrix<double> C = TransformationMatrix<double>::Identity(K, K);

        return SQHamiltonian(ao_basis, S, ScalarSQOneElectronOperator<double>({h}), ScalarSQTwoElectronOperator<double>({g}), C);
    }

    /*
     *  GETTERS
     */

    const ScalarSQOneElectronOperator<Scalar>& get_S() const { return this->S; }
    const TransformationMatrix<Scalar>& get_T_total() const { return this->T_total; }
    const std::shared_ptr<ScalarBasis<GTOShell>>& get_ao_basis() const { return this->ao_basis; }


    /**
     *  @return the dimension of the Hamiltonian, i.e. the number of spinors in which it is expressed
     */
    size_t get_K() const { return this->dimension(); }


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return the dimension of the Hamiltonian, i.e. the number of spinors in which it is expressed
     */
    size_t dimension() const { return this->K; }

    /**
     *  @return the 'core' Hamiltonian, i.e. the total of the one-electron contributions to the Hamiltonian
     */
    const ScalarSQOneElectronOperator<Scalar>& core() const { return this->total_one_op; }

    /**
     *  @return the total of the two-electron contributions to the Hamiltonian
     */
    const ScalarSQTwoElectronOperator<Scalar>& twoElectron() const { return this->total_two_op; }

    /*
     *  PUBLIC METHODS - RELATED TO TRANSFORMATIONS
     */

    /**
     *  @return if the underlying spatial orbital basis of the Hamiltonian is orthonormal
     */
    bool areOrbitalsOrthonormal() const {
        return this->S.parameters().isApprox(SquareMatrix<Scalar>::Identity(this->K, this->K));
    }


    /**
     *  In-place transform the matrix representations of Hamiltonian
     *
     *  @param T    the transformation matrix between the old and the new orbital basis
     *
     *  Furthermore
     *      - the overlap matrix S now gives the overlap matrix in the new molecular orbital basis
     *      - the total transformation matrix T_total is updated to reflect the total transformation between the new molecular orbital basis and the initial atomic orbitals
     */
    void transform(const TransformationMatrix<Scalar>& T) {

        this->S.transform(T);

        this->total_one_op.transform(T);
        this->total_two_op.transform(T);

        this->T_total.transform(T);
    }


    /**
     *  In-place rotate the matrix representations of Hamiltonian
     *
     *  @param U    the unitary rotation matrix between the old and the new orbital basis
     *
     *  Furthermore
     *      - the overlap matrix S now gives the overlap matrix in the new molecular orbital basis
     *      - the total transformation matrix T_total is updated to reflect the total transformation between the new molecular orbital basis and the initial atomic orbitals
     */
    void rotate(const TransformationMatrix<Scalar>& U) {

        this->S.rotate(U);

        this->total_one_op.rotate(U);
        this->total_two_op.rotate(U);

        this->T_total.transform(U);
    }


    /**
     *  In-place rotate the matrix representations of the Hamiltonian using a unitary Jacobi rotation matrix constructed from the Jacobi rotation parameters. Note that this function is only available for real (double) matrix representations
     *
     *  @param jacobi_rotation_parameters       the Jacobi rotation parameters (p, q, angle) that are used to specify a Jacobi rotation: we use the (cos, sin, -sin, cos) definition for the Jacobi rotation matrix
     *
     *  Furthermore
     *      - the overlap matrix S now gives the overlap matrix in the new molecular orbital basis
     *      - the total transformation matrix T_total is updated to reflect the total transformation between the new molecular orbital basis and the initial atomic orbitals
     */
    template<typename Z = Scalar>
    enable_if_t<std::is_same<Z, double>::value> rotate(const JacobiRotationParameters& jacobi_rotation_parameters) {

        this->S.rotate(jacobi_rotation_parameters);
        this->total_one_op.rotate(jacobi_rotation_parameters);
        this->total_two_op.rotate(jacobi_rotation_parameters);

        // Create a Jacobi rotation matrix to transform the coefficient matrix with
        size_t K = this->dimension();  // number of spatial orbitals
        auto J = TransformationMatrix<double>::FromJacobi(jacobi_rotation_parameters, K);
        this->T_total.transform(J);
    }


    /**
     *  Using a random rotation matrix, transform the matrix representations of the Hamiltonian
     *
     *  Note that this method is only available for real matrix representations
     */
    template<typename Z = Scalar>
    enable_if_t<std::is_same<Z, double>::value> randomRotate() {

        // Get a random unitary matrix by diagonalizing a random symmetric matrix
        TransformationMatrix<double> A_random = TransformationMatrix<double>::Random(this->K, this->K);
        TransformationMatrix<double> A_symmetric = A_random + A_random.transpose();
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> unitary_solver (A_symmetric);
        TransformationMatrix<double> U_random = unitary_solver.eigenvectors();

        this->rotate(U_random);
    }


    /**
     *  Transform the SQHamiltonian to the Löwdin basis (i.e. T = S^{-1/2})
     */
    void LowdinOrthonormalize() {

        // The transformation matrix to the Löwdin basis is T = S^{-1/2}
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes (this->S.parameters());
        this->transform(TransformationMatrix<double>(saes.operatorInverseSqrt()));
    }



    /*
     *  PUBLIC METHODS - CALCULATIONS OF VALUES
     */

    /**
     *  @param N_P      the number of electron pairs
     *
     *  @return the Edmiston-Ruedenberg localization index g(i,i,i,i)
     *
     *  Note that this method is only available for real matrix representations
     */
    template<typename Z = Scalar>
    enable_if_t<std::is_same<Z, double>::value, double> calculateEdmistonRuedenbergLocalizationIndex(size_t N_P) const {

        const auto& g_par = this->total_two_op.parameters();

        // TODO: when Eigen releases TensorTrace, use it here
        double localization_index = 0.0;
        for (size_t i = 0; i < N_P; i++) {
            localization_index += g_par(i,i,i,i);
        }

        return localization_index;
    }



    /*
     *  PUBLIC METHODS - CALCULATIONS OF ONE-ELECTRON OPERATORS
     */

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
     *  @param ao_list     indices of the AOs used for the Mulliken populations
     *
     *  @return the Mulliken operator for a set of AOs
     *
     *  Note that this method is only available for real matrix representations
     */
    template<typename Z = Scalar>
    enable_if_t<std::is_same<Z, double>::value, ScalarSQOneElectronOperator<double>> calculateMullikenOperator(const Vectoru& ao_list) const {

        if (!this->get_ao_basis()) {
            throw std::invalid_argument("SQHamiltonian::calculateMullikenOperator(Vectoru): The Hamiltonian has no underlying AO basis, Mulliken analysis is not possible.");
        }

        if (ao_list.size() > this->K) {
            throw std::invalid_argument("SQHamiltonian::calculateMullikenOperator(Vectoru): Too many AOs are selected");
        }

        // Create the partitioning matrix
        SquareMatrix<double> p_a = SquareMatrix<double>::PartitionMatrix(ao_list, this->K);

        ScalarSQOneElectronOperator<Scalar> S_AO = this->S;
        TransformationMatrix<double> T_inverse = T_total.inverse();
        S_AO.transform(T_inverse);

        ScalarSQOneElectronOperator<double> mulliken_matrix ({ (T_total.adjoint() * p_a * S_AO.parameters() * T_total + T_total.adjoint() * S_AO.parameters() * p_a * T_total)/2 });

        return mulliken_matrix;
    }


    /**
     *  @return the effective one-electron integrals
     */
    ScalarSQOneElectronOperator<Scalar> calculateEffectiveOneElectronIntegrals() const {

        return this->core() + this->twoElectron().effectiveOneElectronPartition();
    }



    /*
     *  PUBLIC METHODS - CALCULATIONS OF TWO-ELECTRON OPERATORS
     */

    /**
     *  @param D      the 1-DM (or the response 1-DM for made-variational wave function models)
     *  @param d      the 2-DM (or the response 2-DM for made-variational wave function models)
     *
     *  @return the (generalized) super-Fockian matrix
     */
    SquareRankFourTensor<Scalar> calculateSuperFockianMatrix(const OneRDM<double>& D, const TwoRDM<double>& d) const {

        return this->core().calculateSuperFockianMatrix(D, d)[0].Eigen() + this->twoElectron().calculateSuperFockianMatrix(D, d)[0].Eigen();  // SQHamiltonian are ScalarSQOperators
    }


    /*
     *  PUBLIC METHODS - CONSTRAINTS
     */

    /**
     *  Constrain the Hamiltonian according to the convention: - lambda * constraint
     *
     *  @param one_op   the one-electron operator used as a constraint
     *  @param two_op   the two-electron operator used as a constraint
     *  @param lambda   Lagrangian multiplier for the constraint
     *
     *  @return a copy of the constrained Hamiltonian
     *
     *  Note that this method is only available for real matrix representations
     */
    template<typename Z = Scalar>
    enable_if_t<std::is_same<Z, double>::value, SQHamiltonian<double>> constrain(const ScalarSQOneElectronOperator<double>& one_op, const ScalarSQTwoElectronOperator<double>& two_op, double lambda) const {

        ScalarSQOneElectronOperator<double> h_constrained (this->core() - lambda*one_op);
        ScalarSQTwoElectronOperator<double> g_constrained (this->twoElectron() - lambda*two_op);

        return SQHamiltonian(this->ao_basis, this->S, h_constrained, g_constrained, this->T_total);
    }


    /**
     *  Constrain the Hamiltonian according to the convention: - lambda * constraint
     *
     *  @param one_op   the one-electron operator used as a constraint
     *  @param lambda   Lagrangian multiplier for the constraint
     *
     *  @return a copy of the constrained Hamiltonian
     *
     *  Note that this method is only available for real matrix representations
     */
    template<typename Z = Scalar>
    enable_if_t<std::is_same<Z, double>::value, SQHamiltonian<double>> constrain(const ScalarSQOneElectronOperator<double>& one_op, double lambda) const {

        ScalarSQOneElectronOperator<double> h_constrained (this->core() - lambda*one_op);

        return SQHamiltonian(this->ao_basis, this->S, h_constrained, this->twoElectron(), this->T_total);
    }


    /**
     *  Constrain the Hamiltonian according to the convention: - lambda * constraint
     *
     *  @param two_op   the two-electron operator used as a constraint
     *  @param lambda   Lagrangian multiplier for the constraint
     *
     *  @return a copy of the constrained Hamiltonian
     *
     *  Note that this method is only available for real matrix representations
     */
    template<typename Z = Scalar>
    enable_if_t<std::is_same<Z, double>::value, SQHamiltonian<double>> constrain(const ScalarSQTwoElectronOperator<double>& two_op, double lambda) const {

        ScalarSQTwoElectronOperator<double> g_constrained (this->twoElectron() - lambda*two_op);

        return SQHamiltonian(this->ao_basis, this->S, this->core(), g_constrained, this->T_total);
    }
};


}  // namespace GQCP
