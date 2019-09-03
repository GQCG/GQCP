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
#include "HamiltonianParameters/BaseHamiltonianParameters.hpp"
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
 *  A class for representing Hamiltonian parameters, i.e. the one- and two-electron integrals in the second-quantized expression of the Hamiltonian
 *
 *  This class can be used for restricted calculations, i.e. the alpha and beta integrals are equal
 *
 *  @tparam Scalar      the scalar type
 */
template <typename Scalar>
class HamiltonianParameters : public BaseHamiltonianParameters {
private:
    size_t K;  // the number of spatial orbitals

    ScalarSQOneElectronOperator<Scalar> S;  // overlap

    ScalarSQOneElectronOperator<Scalar> h;  // one-electron interactions (i.e. the core Hamiltonian)
    ScalarSQTwoElectronOperator<Scalar> g;  // two-electron interactions

    TransformationMatrix<Scalar> T_total;  // total transformation matrix between the current (restricted) molecular orbitals and the atomic orbitals


public:

    /*
     *  CONSTRUCTORS
     */
    HamiltonianParameters() = default;

    
    /**
     *  @param ao_basis     the initial AO basis
     *  @param S            the overlap integrals
     *  @param h            the one-electron integrals H_core
     *  @param g            the two-electron integrals
     *  @param C            a transformation matrix between the current molecular orbitals and the atomic orbitals
     *  @param scalar       the scalar interaction term
     */
    HamiltonianParameters(std::shared_ptr<ScalarBasis<GTOShell>> ao_basis, const ScalarSQOneElectronOperator<Scalar>& S, const ScalarSQOneElectronOperator<Scalar>& h, const ScalarSQTwoElectronOperator<Scalar>& g, const TransformationMatrix<Scalar>& C, double scalar=0.0) :
        BaseHamiltonianParameters(std::move(ao_basis), scalar),
        K (S.get_dim()),
        S (S),
        h (h),
        g (g),
        T_total (C)
    {
        // Check if the dimensions of all matrix representations are compatible
        auto error = std::invalid_argument("HamiltonianParameters::HamiltonianParameters(std::shared_ptr<ScalarBasis<GTOShell>>, ScalarSQOneElectronOperator<Scalar>, ScalarSQOneElectronOperator<Scalar>, ScalarSQTwoElectronOperator<Scalar>,TransformationMatrix<Scalar>, double): The dimensions of the operators and coefficient matrix are incompatible.");

        if (this->ao_basis) {  // ao_basis is not nullptr
            if (this->K != this->ao_basis->numberOfBasisFunctions()) {
                throw error;
            }
        }

        if ((h.get_dim() != this->K) || (g.get_dim() != this->K) || (C.cols() != this->K) || (C.rows() != this->K)) {
            throw error;
        }


        if (S.parameters().isZero(1.0e-08)) {
            throw std::invalid_argument("HamiltonianParameters::HamiltonianParameters(std::shared_ptr<ScalarBasis<GTOShell>>, ScalarSQOneElectronOperator<Scalar>, ScalarSQOneElectronOperator<Scalar>, ScalarSQTwoElectronOperator<Scalar>,TransformationMatrix<Scalar>, double): The underlying overlap matrix cannot be a zero matrix.");
        }
    }


    /**
     *  A constructor that transforms the given Hamiltonian parameters with a transformation matrix
     *
     *  @param ham_par      the current Hamiltonian parameters
     *  @param C            the transformation matrix to be applied to the given Hamiltonian parameters
     */
    HamiltonianParameters(const HamiltonianParameters<Scalar>& ham_par, const TransformationMatrix<Scalar>& C) :
        HamiltonianParameters<Scalar>(ham_par.ao_basis, ham_par.S, ham_par.h, ham_par.g, ham_par.T_total, ham_par.scalar)
    {
        // We have now initialized the new Hamiltonian parameters to be a copy of the given Hamiltonian parameters, so now we will transform
        this->transform(C);
    }



    /*
     *  NAMED CONSTRUCTORS
     */

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
    template <typename Z = Scalar>
    static enable_if_t<std::is_same<Z, double>::value, HamiltonianParameters<double>> Molecular(std::shared_ptr<ScalarBasis<GTOShell>> ao_basis, double scalar=0.0) {

        const SingleParticleBasis<Z, GTOShell> sp_basis (*ao_basis);
        const NuclearFramework nuclear_framework (ao_basis->shellSet().nuclei());

        // Calculate the integrals for the molecular Hamiltonian
        const auto S = sp_basis.quantize(Operator::Overlap());
        const auto T = sp_basis.quantize(Operator::Kinetic());
        const auto V = sp_basis.quantize(Operator::NuclearAttraction(nuclear_framework));
        ScalarSQOneElectronOperator<double> H = T + V;

        const auto g = sp_basis.quantize(Operator::Coulomb());


        // Construct the initial transformation matrix: the identity matrix
        auto nbf = ao_basis->numberOfBasisFunctions();
        TransformationMatrix<double> T_total = TransformationMatrix<double>::Identity(nbf, nbf);

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

        auto ao_basis = std::make_shared<ScalarBasis<GTOShell>>(molecule, basisset);

        const double internuclear_repulsion_energy = Operator::NuclearRepulsion(molecule).value();
        return HamiltonianParameters::Molecular(ao_basis, internuclear_repulsion_energy);
    }


    /**
     *  @param K        the number of orbitals
     *
     *  @return a set of random Hamiltonian parameters with values uniformly distributed between [-1,1]
     *
     *  Note that this named constructor is only available for real representations
     */
    template <typename Z = Scalar>
    static enable_if_t<std::is_same<Z, double>::value, HamiltonianParameters<double>> Random(size_t K) {


        ScalarSQOneElectronOperator<double> S ({ChemicalMatrix<double>::Identity(K, K)});  // the underlying orbital basis can be chosen as orthonormal, since the form of the underlying orbitals doesn't really matter
        TransformationMatrix<double> C = TransformationMatrix<double>::Identity(K, K);  // the transformation matrix here doesn't really mean anything, because it doesn't link to any AO basis

        ScalarSQOneElectronOperator<double> H ({ChemicalMatrix<double>::Random(K, K)});  // uniformly distributed between [-1,1]


        // Unfortunately, the Tensor module provides uniform random distributions between [0, 1]
        ChemicalRankFourTensor<double> g (K);
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


        // Get a random scalar
        std::random_device random_device;  // used to seed PRNG
        std::mt19937 random_generator (random_device());
        std::uniform_real_distribution<double> double_distribution (-1.0, 1.0);
        double scalar = double_distribution(random_generator);

        return HamiltonianParameters<double>(ao_basis, S, H, ScalarSQTwoElectronOperator<double>({g}), C, scalar);
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
            throw std::invalid_argument("HamiltonianParameters::ReadFCIDUMP(std::string): The .FCIDUMP-file is invalid: could not read a number of orbitals.");
        }


        double scalar = 0.0;
        ChemicalMatrix<double> h_core = ChemicalMatrix<double>::Zero(K, K);
        ChemicalRankFourTensor<double> g (K);
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


        // Make the ingredients to construct HamiltonianParameters
        std::shared_ptr<ScalarBasis<GTOShell>> ao_basis;  // nullptr
        ScalarSQOneElectronOperator<Scalar> S ({ChemicalMatrix<double>::Identity(K, K)});
        TransformationMatrix<double> C = TransformationMatrix<double>::Identity(K, K);

        return HamiltonianParameters(ao_basis, S, ScalarSQOneElectronOperator<Scalar>({h_core}), ScalarSQTwoElectronOperator<Scalar>({g}), C, scalar);
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

        const size_t K = H.numberOfLatticeSites();

        ChemicalMatrix<double> h = ChemicalMatrix<double>::Zero(K, K);
        ChemicalRankFourTensor<double> g (K);
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


        // Make the ingredients to construct HamiltonianParameters
        std::shared_ptr<ScalarBasis<GTOShell>> ao_basis;  // nullptr
        ScalarSQOneElectronOperator<double> S ({ChemicalMatrix<double>::Identity(K, K)});
        TransformationMatrix<double> C = TransformationMatrix<double>::Identity(K, K);

        return HamiltonianParameters(ao_basis, S, ScalarSQOneElectronOperator<double>({h}), ScalarSQTwoElectronOperator<double>({g}), C);  // no scalar term
    }



    /*
     *  DESTRUCTOR
     */

    ~HamiltonianParameters() override = default;


    /*
     *  GETTERS
     */

    const ScalarSQOneElectronOperator<Scalar>& get_S() const { return this->S; }
    const ScalarSQOneElectronOperator<Scalar>& get_h() const { return this->h; }
    const ScalarSQTwoElectronOperator<Scalar>& get_g() const { return this->g; }
    const TransformationMatrix<Scalar>& get_T_total() const { return this->T_total; }
    size_t get_K() const { return this->K; }

    size_t dimension() const { return this->K; }


    /*
     *  PUBLIC METHODS - RELATED TO TRANSFORMATIONS
     */

    /**
     *  @return if the underlying spatial orbital basis of the Hamiltonian parameters is orthonormal
     */
    bool areOrbitalsOrthonormal() const {
        return this->S.parameters().isApprox(SquareMatrix<Scalar>::Identity(this->K, this->K));
    }


    /**
     *  In-place transform the matrix representations of Hamiltonian parameters
     *
     *  @param T    the transformation matrix between the old and the new orbital basis
     *
     *  Furthermore
     *      - the overlap matrix S now gives the overlap matrix in the new molecular orbital basis
     *      - the total transformation matrix T_total is updated to reflect the total transformation between the new molecular orbital basis and the initial atomic orbitals
     */
    void transform(const TransformationMatrix<Scalar>& T) {

        this->S.transform(T);

        this->h.transform(T);
        this->g.transform(T);

        this->T_total.transform(T);
    }


    /**
     *  In-place rotate the matrix representations of Hamiltonian parameters
     *
     *  @param U    the unitary rotation matrix between the old and the new orbital basis
     *
     *  Furthermore
     *      - the overlap matrix S now gives the overlap matrix in the new molecular orbital basis
     *      - the total transformation matrix T_total is updated to reflect the total transformation between the new molecular orbital basis and the initial atomic orbitals
     */
    void rotate(const TransformationMatrix<Scalar>& U) {

        this->S.rotate(U);

        this->h.rotate(U);
        this->g.rotate(U);

        this->T_total.transform(U);
    }


    /**
     *  In-place rotate the matrix representations of the Hamiltonian parameters using a unitary Jacobi rotation matrix constructed from the Jacobi rotation parameters. Note that this function is only available for real (double) matrix representations
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
        this->h.rotate(jacobi_rotation_parameters);
        this->g.rotate(jacobi_rotation_parameters);

        // Create a Jacobi rotation matrix to transform the coefficient matrix with
        size_t K = this->h.get_dim();  // number of spatial orbitals
        auto J = TransformationMatrix<double>::FromJacobi(jacobi_rotation_parameters, K);
        this->T_total.transform(J);
    }


    /**
     *  Using a random rotation matrix, transform the matrix representations of the Hamiltonian parameters
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
     *  Transform the HamiltonianParameters to the Löwdin basis (i.e. T = S^{-1/2})
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

        const auto& g_par = this->g.parameters();

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

        return this->h.calculateFockianMatrix(D, d)[0] + this->g.calculateFockianMatrix(D, d)[0];  // HamiltonianParameters are a combination of scalar parameters
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
            throw std::invalid_argument("HamiltonianParameters::calculateMullikenOperator(Vectoru): The Hamiltonian parameters have no underlying AO basis, Mulliken analysis is not possible.");
        }

        if (ao_list.size() > this->K) {
            throw std::invalid_argument("HamiltonianParameters::calculateMullikenOperator(Vectoru): Too many AOs are selected");
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

        return this->h + this->g.effectiveOneElectronPartition();
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

        return this->h.calculateSuperFockianMatrix(D, d)[0].Eigen() + this->g.calculateSuperFockianMatrix(D, d)[0].Eigen();  // HamiltonianParameters are ScalarSQOperators
    }


    /*
     *  PUBLIC METHODS - CONSTRAINTS
     */

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
    enable_if_t<std::is_same<Z, double>::value, HamiltonianParameters<double>> constrain(const ScalarSQOneElectronOperator<double>& one_op, const ScalarSQTwoElectronOperator<double>& two_op, double lambda) const {

        ScalarSQOneElectronOperator<double> h_constrained (this->h - lambda*one_op);
        ScalarSQTwoElectronOperator<double> g_constrained (this->g - lambda*two_op);

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
    enable_if_t<std::is_same<Z, double>::value, HamiltonianParameters<double>> constrain(const ScalarSQOneElectronOperator<double>& one_op, double lambda) const {

        ScalarSQOneElectronOperator<double> h_constrained (this->h - lambda*one_op);

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
    enable_if_t<std::is_same<Z, double>::value, HamiltonianParameters<double>> constrain(const ScalarSQTwoElectronOperator<double>& two_op, double lambda) const {

        ScalarSQTwoElectronOperator<double> g_constrained (this->g - lambda*two_op);

        return HamiltonianParameters(this->ao_basis, this->S, this->h, g_constrained, this->T_total);
    }
};


}  // namespace GQCP
