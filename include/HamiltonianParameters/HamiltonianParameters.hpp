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
#ifndef GQCP_HAMILTONIANPARAMETERS_HPP
#define GQCP_HAMILTONIANPARAMETERS_HPP

#include "HamiltonianParameters/BaseHamiltonianParameters.hpp"
#include "HoppingMatrix.hpp"
#include "JacobiRotationParameters.hpp"
#include "LibintCommunicator.hpp"
#include "Molecule.hpp"
#include "Operator/OneElectronOperator.hpp"
#include "Operator/TwoElectronOperator.hpp"
#include "RDM/TwoRDM.hpp"
#include "RDM/OneRDM.hpp"
#include "typedefs.hpp"


namespace GQCP {


/**
 *  A class for representing Hamiltonian parameters, i.e. the one- and two-electron integrals in the second-quantized expression of the Hamiltonian
 *
 *  This class can be used for restricted calculations, i.e. the alpha and beta integrals are equal
 *
 *  @tparam Scalar      the scalar type
 */
template<typename Scalar>
class HamiltonianParameters : public BaseHamiltonianParameters, public Operator<Scalar> {
private:
    size_t K;  // the number of spatial orbitals

    OneElectronOperator<Scalar> S;  // overlap

    OneElectronOperator<Scalar> h;  // one-electron interactions (i.e. the core Hamiltonian)
    TwoElectronOperator<Scalar> g;  // two-electron interactions

    SquareMatrix<Scalar> T_total;  // total transformation matrix between the current (restricted) molecular orbitals and the atomic orbitals


public:
    // CONSTRUCTORS
    /**
     *  @param ao_basis     the initial AO basis
     *  @param S            the overlap integrals
     *  @param h            the one-electron integrals H_core
     *  @param g            the two-electron integrals
     *  @param C            a transformation matrix between the current molecular orbitals and the atomic orbitals
     *  @param scalar       the scalar interaction term
     */
    HamiltonianParameters(std::shared_ptr<AOBasis> ao_basis, const OneElectronOperator<Scalar>& S, const OneElectronOperator<Scalar>& h, const TwoElectronOperator<Scalar>& g, const SquareMatrix<Scalar>& C, double scalar=0.0) :
        BaseHamiltonianParameters(std::move(ao_basis), scalar),
        K (S.get_dim()),
        S (S),
        h (h),
        g (g),
        T_total (C)
    {
        // Check if the dimensions of all matrix representations are compatible
        auto error = std::invalid_argument("The dimensions of the operators and coefficient matrix are incompatible.");

        if (this->ao_basis) {  // ao_basis is not nullptr
            if (this->K != this->ao_basis->get_number_of_basis_functions()) {
                throw error;
            }
        }

        if ((h.get_dim() != this->K) || (g.get_dim() != this->K) || (C.cols() != this->K) || (C.rows() != this->K)) {
            throw error;
        }


        if (S.isZero(1.0e-08)) {
            throw std::invalid_argument("The underlying overlap matrix cannot be a zero matrix.");
        }
    }


    /**
     *  A constructor that transforms the given Hamiltonian parameters with a transformation matrix
     *
     *  @param ham_par      the current Hamiltonian parameters
     *  @param C            the transformation matrix to be applied to the given Hamiltonian parameters
     */
    HamiltonianParameters(const HamiltonianParameters<Scalar>& ham_par, const SquareMatrix<Scalar>& C) :
        HamiltonianParameters<Scalar>(ham_par.ao_basis, ham_par.S, ham_par.h, ham_par.g, ham_par.T_total, ham_par.scalar)
    {
        // We have now initialized the new Hamiltonian parameters to be a copy of the given Hamiltonian parameters, so now we will transform
        this->transform(C);
    }


    // NAMED CONSTRUCTORS
    /**
     *  Construct the molecular Hamiltonian parameters in an AO basis
     *
     *  @param ao_basis     the AO basis in which the Hamiltonian parameters should be expressed
     *  @param scalar       the scalar energy term (usually the internuclear repulsion energy)
     *
     *  @return Hamiltonian parameters corresponding to the molecular Hamiltonian in an AO basis. The molecular Hamiltonian has
     *      - scalar contributions:
     *          - internuclear repulsion
     *      - one-electron contributions:
     *          - kinetic
     *          - nuclear attraction
     *      - two-electron contributions:
     *          - Coulomb repulsion
     *
     *  Note that this named constructor is only available for real matrix representations
     */
    template<typename Z = Scalar>
    static enable_if_t<std::is_same<Z, double>::value, HamiltonianParameters<double>> Molecular(std::shared_ptr<AOBasis> ao_basis, double scalar=0.0) {

        // Calculate the integrals for the molecular Hamiltonian
        auto S = LibintCommunicator::get().calculateOverlapIntegrals(*ao_basis);
        auto T = LibintCommunicator::get().calculateKineticIntegrals(*ao_basis);
        auto V = LibintCommunicator::get().calculateNuclearIntegrals(*ao_basis);
        auto H = OneElectronOperator<double>(T + V);

        auto g = LibintCommunicator::get().calculateCoulombRepulsionIntegrals(*ao_basis);


        // Construct the initial transformation matrix: the identity matrix
        auto nbf = ao_basis->get_number_of_basis_functions();
        auto T_total = SquareMatrix<double>(Matrix<double>::Identity(nbf, nbf));

        return HamiltonianParameters(ao_basis, S, H, g, T_total, scalar);
    }


    /**
     *  Construct the molecular Hamiltonian parameters in an AO basis
     *
     *  @param molecule     the molecule for which the Hamiltonian parameters should be calculated
     *  @param basisset     the name of the basisset corresponding to the AO basis
     *
     *  @return Hamiltonian parameters corresponding to the molecular Hamiltonian in an AO basis. The molecular Hamiltonian has
     *      - scalar contributions:
     *          - internuclear repulsion
     *      - one-electron contributions:
     *          - kinetic
     *          - nuclear attraction
     *      - two-electron contributions:
     *          - Coulomb repulsion
     *
     *  Note that this named constructor is only available for real matrix representations
     */
    template<typename Z = Scalar>
    static enable_if_t<std::is_same<Z, double>::value, HamiltonianParameters<double>> Molecular(const Molecule& molecule, const std::string& basisset) {

        auto ao_basis = std::make_shared<AOBasis>(molecule, basisset);

        return HamiltonianParameters::Molecular(ao_basis, molecule.calculateInternuclearRepulsionEnergy());
    }

    /**
     *  @param K        the number of orbitals
     *
     *  @return a set of random Hamiltonian parameters with values uniformly distributed between [-1,1]
     *
     *  Note that this named constructor is only available for real representations
     */
    template<typename Z = Scalar>
    static enable_if_t<std::is_same<Z, double>::value, HamiltonianParameters<double>> Random(size_t K) {
        OneElectronOperator<double> S (Matrix<double>::Identity(K, K));  // the underlying orbital basis can be chosen as orthonormal, since the form of the underlying orbitals doesn't really matter
        Matrix<double> C (Matrix<double>::Identity(K, K));  // the transformation matrix C here doesn't really mean anything, because it doesn't link to any AO basis

        OneElectronOperator<double> H (Matrix<double>::Random(K, K));  // uniformly distributed between [-1,1]


        // Unfortunately, the Tensor module provides uniform random distributions between [0, 1]
        FourIndexTensor<double> g_tensor (K, K, K, K);
        g_tensor.setRandom();

        // Move the distribution from [0, 1] -> [-1, 1]
        for (size_t i = 0; i < K; i++) {
            for (size_t j = 0; j < K; j++) {
                for (size_t k = 0; k < K; k++) {
                    for (size_t l = 0; l < K; l++) {
                        g_tensor(i,j,k,l) = 2*g_tensor(i,j,k,l) - 1;  // scale from [0, 1] -> [0, 2] -> [-1, 1]
                    }
                }
            }
        }
        TwoElectronOperator<double> g (g_tensor);

        std::shared_ptr<AOBasis> ao_basis;  // nullptr because it doesn't make sense to set an AOBasis


        // Get a random scalar
        std::random_device random_device;  // used to seed PRNG
        std::mt19937 random_generator (random_device());
        std::uniform_real_distribution<double> double_distribution (-1.0, 1.0);
        double scalar = double_distribution(random_generator);

        return HamiltonianParameters<double>(ao_basis, S, H, g, C, scalar);
    }

    /**
     *  @param fcidump_file     the name of the FCIDUMP file
     *
     *  @return Hamiltonian parameters corresponding to the contents of an FCIDUMP file
     *
     *  Note that this named constructor is only available for real matrix representations
     */
    template<typename Z = Scalar>
    static enable_if_t<std::is_same<Z, double>::value, HamiltonianParameters<double>> ReadFCIDUMP(const std::string& fcidump_file) {

        // Find the extension of the given path (https://stackoverflow.com/a/51992)
        std::string extension;
        std::string::size_type idx = fcidump_file.rfind('.');

        if (idx != std::string::npos) {
            extension = fcidump_file.substr(idx+1);
        } else {
            throw std::runtime_error("I did not find an extension in your given path.");
        }

        if (!(extension == "FCIDUMP")) {
            throw std::runtime_error("You did not provide a .FCIDUMP file name");
        }

        // If the xyz_filename isn't properly converted into an input file stream, we assume the user supplied a wrong file
        std::ifstream input_file_stream (fcidump_file);

        if (!input_file_stream.good()) {
            throw std::runtime_error("The provided FCIDUMP file is illegible. Maybe you specified a wrong path?");
        }



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
            throw std::invalid_argument("The .FCIDUMP-file is invalid: could not read a number of orbitals.");
        }


        double scalar = 0.0;
        Eigen::MatrixXd h_SO = Eigen::MatrixXd::Zero(K, K);
        Eigen::Tensor<double, 4> g_SO (K, K, K, K);
        g_SO.setZero();

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

            //  Internuclear repulsion energy
            if ((i == 0) && (j == 0) && (a == 0) && (b == 0)) {
                scalar = x;
            }

            //  Single-particle eigenvalues (skipped)
            else if ((a == 0) && (j == 0) && (b == 0)) {}

            //  One-electron integrals (h_core)
            else if ((j == 0) && (b == 0)) {
                size_t p = i - 1;
                size_t q = a - 1;
                h_SO(p,q) = x;

                // Apply the permutational symmetry for real orbitals
                h_SO(q,p) = x;
            }

            //  Two-electron integrals are given in CHEMIST'S NOTATION, so just copy them over
            else if ((i > 0) && (a > 0) && (j > 0) && (b > 0)) {
                size_t p = i - 1;
                size_t q = a - 1;
                size_t r = j - 1;
                size_t s = b - 1;
                g_SO(p,q,r,s) = x;

                // Apply the permutational symmetries for real orbitals
                g_SO(p,q,s,r) = x;
                g_SO(q,p,r,s) = x;
                g_SO(q,p,s,r) = x;

                g_SO(r,s,p,q) = x;
                g_SO(s,r,p,q) = x;
                g_SO(r,s,q,p) = x;
                g_SO(s,r,q,p) = x;
            }
        }  // while loop


        // Make the ingredients to construct HamiltonianParameters
        std::shared_ptr<AOBasis> ao_basis;  // nullptr
        OneElectronOperator<Scalar> S (Eigen::MatrixXd::Identity(K, K));
        OneElectronOperator<Scalar> H_core (h_SO);
        TwoElectronOperator<Scalar> G (g_SO);
        Eigen::MatrixXd C = Eigen::MatrixXd::Identity(K, K);

        return HamiltonianParameters(ao_basis, S, H_core, G, C, scalar);
    }

    /**
     *  @param H      a Hubbard hopping matrix
     *
     *  @return Hubbard Hamiltonian parameters generated from the Hubbard hopping matrix
     *
     *  Note that this named constructor is only available for real matrix representations
     */
    template<typename Z = Scalar>
    static enable_if_t<std::is_same<Z, double>::value, HamiltonianParameters<double>> Hubbard(const HoppingMatrix& H) {

        size_t K = H.numberOfLatticeSites();

        Eigen::MatrixXd h_SO = Eigen::MatrixXd::Zero(K, K);
        Eigen::Tensor<double, 4> g_SO (K, K, K, K);
        g_SO.setZero();


        for (size_t i = 0; i < K; i++) {
            for (size_t j = i; j < K; j++) {
                if (i == j) {
                    g_SO(i,i,i,i) = H(i,i);
                } else {
                    h_SO(i,j) = H(i,j);
                    h_SO(j,i) = H(j,i);
                }
            }
        }


        // Make the ingredients to construct HamiltonianParameters
        std::shared_ptr<AOBasis> ao_basis;  // nullptr
        OneElectronOperator<Scalar> S (Eigen::MatrixXd::Identity(K, K));
        OneElectronOperator<Scalar> H_core (h_SO);
        TwoElectronOperator<Scalar> G (g_SO);
        Eigen::MatrixXd C = Eigen::MatrixXd::Identity(K, K);

        return HamiltonianParameters(ao_basis, S, H_core, G, C);  // no scalar term
    }


    // DESTRUCTORS
    ~HamiltonianParameters() override = default;

    
    // GETTERS
    const OneElectronOperator<Scalar>& get_S() const { return this->S; }
    const OneElectronOperator<Scalar>& get_h() const { return this->h; }
    const TwoElectronOperator<Scalar>& get_g() const { return this->g; }
    const SquareMatrix<Scalar>& get_T_total() const { return this->T_total; }
    size_t get_K() const { return this->K; }


    // PUBLIC METHODS
    /**
     *  @return if the underlying spatial orbital basis of the Hamiltonian parameters is orthonormal
     */
    bool areOrbitalsOrthonormal() const {
        return this->S.isApprox(Matrix<Scalar>::Identity(this->K, this->K));
    }


    // OVERRIDDEN PUBLIC METHODS
    /**
     *  In-place transform the matrix representations of Hamiltonian parameters
     *
     *  @param T    the transformation matrix between the old and the new orbital basis, it is used as
     *      b' = b T ,
     *   in which the basis functions are collected as elements of a row vector b
     *
     *  Furthermore
     *      - the overlap matrix S now gives the overlap matrix in the new molecular orbital basis
     *      - the total transformation matrix T_total is updated to reflect the total transformation between the new molecular orbital basis and the initial atomic orbitals
     */
    void transform(const SquareMatrix<Scalar>& T) override {

        this->S.transform(T);

        this->h.transform(T);
        this->g.transform(T);

        this->T_total = SquareMatrix<Scalar>(this->T_total * T);  // use the correct transformation formula for subsequent transformations
    }


    /**
     *  In-place rotate the matrix representations of the Hamiltonian parameters using a unitary Jacobi rotation matrix constructed from the Jacobi rotation parameters. Note that this function is only available for real (double) matrix representations
     *
     *  @param jacobi_rotation_parameters       the Jacobi rotation parameters (p, q, angle) that are used to specify a Jacobi rotation: we use the (cos, sin, -sin, cos) definition for the Jacobi rotation matrix. See transform() for how the transformation matrix between the two bases should be represented
     *
     *  Furthermore
     *      - the overlap matrix S now gives the overlap matrix in the new molecular orbital basis
     *      - the total transformation matrix T_total is updated to reflect the total transformation between the new molecular orbital basis and the initial atomic orbitals
     */
    template<typename Z = Scalar>
    enable_if_t<std::is_same<Z, double>::value> rotate(const JacobiRotationParameters& jacobi_rotation_parameters) {

        this->S.rotate(jacobi_rotation_parameters);
        this->h.rotate(jacobi_rotation_parameters);
        this->g.rotate(jacobi_rotation_parameters);


        // Create a Jacobi rotation matrix to transform the coefficient matrix with
        size_t K = this->h.get_dim();  // number of spatial orbitals
        auto J = jacobiRotationMatrix(jacobi_rotation_parameters, K);
        this->T_total = SquareMatrix<Scalar>(this->T_total * J);
    }


    // OTHER TRANSFORMATION FORMULAS
    /**
     *  Using a random rotation matrix, transform the matrix representations of the Hamiltonian parameters
     *
     *  Note that this method is only available for real matrix representations
     */
    template<typename Z = Scalar>
    enable_if_t<std::is_same<Z, double>::value> randomRotate() {

        // Get a random unitary matrix by diagonalizing a random symmetric matrix
        Eigen::MatrixXd A_random = Eigen::MatrixXd::Random(this->K, this->K);
        Eigen::MatrixXd A_symmetric = A_random + A_random.transpose();
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> unitary_solver (A_symmetric);
        Eigen::MatrixXd U_random = unitary_solver.eigenvectors();

        this->rotate(U_random);
    }

    /**
     *  Transform the HamiltonianParameters to the Löwdin basis (i.e. T = S^{-1/2})
     */
    void LowdinOrthonormalize() {

        // The transformation matrix to the Löwdin basis is T = S^{-1/2}
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes (this->S);  // can we use this->S?
        this->transform(saes.operatorInverseSqrt());
    }


    // PUBLIC METHODS - CALCULATIONS OF VALUES
    /**
     *  @param N_P      the number of electron pairs
     *
     *  @return the Edmiston-Ruedenberg localization index g(i,i,i,i)
     *
     *  Note that this method is only available for real matrix representations
     */
    template<typename Z = Scalar>
    enable_if_t<std::is_same<Z, double>::value, double> calculateEdmistonRuedenbergLocalizationIndex(size_t N_P) const {

        double localization_index = 0.0;

        // TODO: when Eigen releases TensorTrace, use it here
        for (size_t i = 0; i < N_P; i++) {
            localization_index += this->g(i,i,i,i);
        }

        return localization_index;
    }


    // PUBLIC METHODS - CALCULATIONS OF ONE-ELECTRON OPERATORS
    /**
     *  @param D      the 1-RDM
     *  @param d      the 2-RDM
     *
     *  @return the generalized Fock matrix
     */
    OneElectronOperator<Scalar> calculateGeneralizedFockMatrix(const OneRDM& D, const TwoRDM& d) const {

        // Check if dimensions are compatible
        if (D.get_dim() != this->K) {
            throw std::invalid_argument("The 1-RDM is not compatible with the HamiltonianParameters.");
        }

        if (d.get_dim() != this->K) {
            throw std::invalid_argument("The 2-RDM is not compatible with the HamiltonianParameters.");
        }


        // A KISS implementation of the calculation of the generalized Fock matrix F
        Eigen::MatrixXd F = Eigen::MatrixXd::Zero(this->K, this->K);
        for (size_t p = 0; p < this->K; p++) {
            for (size_t q = 0; q < this->K; q++) {

                // One-electron part
                for (size_t r = 0; r < this->K; r++) {
                    F(p,q) += h(q,r) * D(p,r);
                }

                // Two-electron part
                for (size_t r = 0; r < this->K; r++) {
                    for (size_t s = 0; s < this->K; s++) {
                        for (size_t t = 0; t < this->K; t++) {
                            F(p,q) += g(q,r,s,t) * d(p,r,s,t);
                        }
                    }
                }  // two-electron part

            }
        }  // F elements loop


        return OneElectronOperator<Scalar>(F);
    }

    /**
     *  @param ao_list     indices of the AOs used for the Mulliken populations
     *
     *  @return the Mulliken operator for a set of AOs
     *
     *  Note that this method is only available for real matrix representations
     */
    template<typename Z = Scalar>
    enable_if_t<std::is_same<Z, double>::value, OneElectronOperator<double>> calculateMullikenOperator(const Vectoru& ao_list) const {

        if (!this->get_ao_basis()) {
            throw std::invalid_argument("The Hamiltonian parameters have no underlying AO basis, Mulliken analysis is not possible.");
        }

        if (ao_list.size() > this->K) {
            throw std::invalid_argument("Too many AOs are selected");
        }

        // Create the partitioning matrix (diagonal matrix with values set to 1 of selected AOs
        Eigen::MatrixXd p_a = Eigen::MatrixXd::Zero(this->K, this->K);

        for (size_t index : ao_list) {
            if (index >= this->K) {
                throw std::invalid_argument("AO index is too large");
            }

            p_a(index, index) = 1;
        }

        OneElectronOperator<Scalar> S_AO = this->S;
        auto T_inverse = SquareMatrix<double>(T_total.inverse());
        S_AO.transform(T_inverse);

        Eigen::MatrixXd mulliken_matrix = (T_total.adjoint() * p_a * S_AO * T_total + T_total.adjoint() * S_AO * p_a * T_total)/2 ;

        return OneElectronOperator<double>(mulliken_matrix);
    }

    /**
     *  @return the effective one-electron integrals
     */
    OneElectronOperator<Scalar> calculateEffectiveOneElectronIntegrals() const {

        auto k = this->h;

        for (size_t p = 0; p < this->K; p++) {
            for (size_t q = 0; q < this->K; q++) {
                for (size_t r = 0; r < this->K; r++) {
                    k(p,q) -= 0.5 * this->g(p,r,r,q);
                }
            }
        }

        return k;
    }


    // PUBLIC METHODS - CALCULATIONS OF TWO-ELECTRON OPERATORS
    /**
     *  @param D      the 1-RDM
     *  @param d      the 2-RDM
     *
     *  @return the super-generalized Fock matrix
     */
    TwoElectronOperator<Scalar> calculateSuperGeneralizedFockMatrix(const OneRDM& D, const TwoRDM& d) const {

        // Check if dimensions are compatible
        if (D.get_dim() != this->K) {
            throw std::invalid_argument("The 1-RDM is not compatible with the HamiltonianParameters.");
        }

        if (d.get_dim() != this->K) {
            throw std::invalid_argument("The 2-RDM is not compatible with the HamiltonianParameters.");
        }


        // We have to calculate the generalized Fock matrix F first
        OneElectronOperator<Scalar> F = this->calculateGeneralizedFockMatrix(D, d);

        // A KISS implementation of the calculation of the super generalized Fock matrix W
        Eigen::Tensor<double, 4> W (this->K, this->K, this->K, this->K);
        W.setZero();
        for (size_t p = 0; p < this->K; p++) {
            for (size_t q = 0; q < this->K; q++) {
                for (size_t r = 0; r < this->K; r++) {
                    for (size_t s = 0; s < this->K; s++) {

                        // Generalized Fock matrix part
                        if (r == q) {
                            W(p,q,r,s) += F(p,s);
                        }

                        // One-electron part
                        W(p,q,r,s) -= this->h(s,p) * D(r,q);

                        // Two-electron part
                        for (size_t t = 0; t < this->K; t++) {
                            for (size_t u = 0; u < this->K; u++) {
                                W(p,q,r,s) += this->g(s,t,q,u) * d(r,t,p,u) - this->g(s,t,u,p) * d(r,t,u,q) - this->g(s,p,t,u) * d(r,q,t,u);
                            }
                        }  // two-electron part
                    }
                }
            }
        }  // W elements loop


        return TwoElectronOperator<Scalar>(W);
    }


    // PUBLIC METHODS - CONSTRAINTS
    /**
     *  Constrain the Hamiltonian parameters according to the convention: - lambda * constraint
     *
     *  @param one_op   the one-electron operator used as a constraint
     *  @param two_op   the two-electron operator used as a constraint
     *  @param lambda   Lagrangian multiplier for the constraint
     *
     *  @return a copy of the constrained Hamiltonian parameters
     *
     *  Note that this method is only available for real matrix representations
     */
    template<typename Z = Scalar>
    enable_if_t<std::is_same<Z, double>::value, HamiltonianParameters<double>> constrain(const OneElectronOperator<double>& one_op, const TwoElectronOperator<double>& two_op, double lambda) const {

        OneElectronOperator<double> h_constrained (this->h - lambda*one_op);
        TwoElectronOperator<double> g_constrained (this->g - lambda*two_op);

        return HamiltonianParameters(this->ao_basis, this->S, h_constrained, g_constrained, this->T_total);
    }

    /**
     *  Constrain the Hamiltonian parameters according to the convention: - lambda * constraint
     *
     *  @param one_op   the one-electron operator used as a constraint
     *  @param lambda   Lagrangian multiplier for the constraint
     *
     *  @return a copy of the constrained Hamiltonian parameters
     *
     *  Note that this method is only available for real matrix representations
     */
    template<typename Z = Scalar>
    enable_if_t<std::is_same<Z, double>::value, HamiltonianParameters<double>> constrain(const OneElectronOperator<double>& one_op, double lambda) const {

        OneElectronOperator<double> h_constrained (this->h - lambda*one_op);

        return HamiltonianParameters(this->ao_basis, this->S, h_constrained, this->g, this->T_total);
    }

    /**
     *  Constrain the Hamiltonian parameters according to the convention: - lambda * constraint
     *
     *  @param two_op   the two-electron operator used as a constraint
     *  @param lambda   Lagrangian multiplier for the constraint
     *
     *  @return a copy of the constrained Hamiltonian parameters
     *
     *  Note that this method is only available for real matrix representations
     */
    template<typename Z = Scalar>
    enable_if_t<std::is_same<Z, double>::value, HamiltonianParameters<double>> constrain(const TwoElectronOperator<double>& two_op, double lambda) const {

        TwoElectronOperator<double> g_constrained (this->g - lambda*two_op);

        return HamiltonianParameters(this->ao_basis, this->S, this->h, g_constrained, this->T_total);
    }
};


}  // namespace GQCP


#endif  // GQCP_HAMILTONIANPARAMETERS_HPP
