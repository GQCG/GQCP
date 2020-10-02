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
#include "Basis/SpinorBasis/GSpinorBasis.hpp"
#include "Basis/SpinorBasis/OrbitalSpace.hpp"
#include "Basis/SpinorBasis/RSpinorBasis.hpp"
#include "Basis/Transformations/JacobiRotationParameters.hpp"
#include "Basis/Transformations/TransformationMatrix.hpp"
#include "DensityMatrix/OneDM.hpp"
#include "DensityMatrix/TwoDM.hpp"
#include "Molecule/Molecule.hpp"
#include "Operator/FirstQuantized/NuclearRepulsionOperator.hpp"
#include "Operator/FirstQuantized/OverlapOperator.hpp"
#include "Operator/SecondQuantized/ModelHamiltonian/HubbardHamiltonian.hpp"
#include "Operator/SecondQuantized/SQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/SQTwoElectronOperator.hpp"
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

        const auto dim = one_ops[0].numberOfOrbitals();
        for (const auto& one_op : this->one_ops) {
            if (one_op.numberOfOrbitals() != dim) {
                throw dimension_error;
            }
        }

        for (const auto& two_op : this->two_ops) {
            if (two_op.numberOfOrbitals() != dim) {
                throw dimension_error;
            }
        }


        // Calculate the total one-electron operator
        SquareMatrix<Scalar> total_one_op_par {dim};
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
        SQHamiltonian(std::vector<ScalarSQOneElectronOperator<Scalar>> {h}, std::vector<ScalarSQTwoElectronOperator<Scalar>> {g}) {}


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
     *  Construct the molecular Hamiltonian in a given restricted spin-orbital basis.
     *
     *  @param r_spinor_basis       the spinor basis in which the Hamiltonian should be expressed
     *  @param molecule             the molecule on which the single particle is based
     *
     *  @return a second-quantized molecular Hamiltonian. The molecular Hamiltonian has
     *      - one-electron contributions:
     *          - kinetic
     *          - nuclear attraction
     *      - two-electron contributions:
     *          - Coulomb repulsion
     *
     *  @note This named constructor is only available for real matrix representations.
     */
    template <typename Z = Scalar>
    static enable_if_t<std::is_same<Z, double>::value, SQHamiltonian<double>> Molecular(const RSpinorBasis<Z, GTOShell>& r_spinor_basis, const Molecule& molecule) {

        // Calculate the integrals for the molecular Hamiltonian
        const auto T = r_spinor_basis.quantize(Operator::Kinetic());
        const auto V = r_spinor_basis.quantize(Operator::NuclearAttraction(molecule));
        ScalarSQOneElectronOperator<double> H = T + V;

        const auto g = r_spinor_basis.quantize(Operator::Coulomb());

        return SQHamiltonian(H, g);
    }


    /**
     *  Construct the molecular Hamiltonian in a given (general) spinor basis.
     *
     *  @param g_spinor_basis           the (general) spinor basis in which the Hamiltonian should be expressed
     *  @param molecule                 the molecule on which the single particle is based
     *
     *  @return a second-quantized molecular Hamiltonian. The molecular Hamiltonian has
     *      - one-electron contributions:
     *          - kinetic
     *          - nuclear attraction
     *      - two-electron contributions:
     *          - Coulomb repulsion
     *
     *  @note This named constructor is only available for real matrix representations.
     */
    template <typename Z = Scalar>
    static enable_if_t<std::is_same<Z, double>::value, SQHamiltonian<double>> Molecular(const GSpinorBasis<Z, GTOShell>& g_spinor_basis, const Molecule& molecule) {

        // Calculate the integrals for the molecular Hamiltonian
        const auto T = g_spinor_basis.quantize(Operator::Kinetic());
        const auto V = g_spinor_basis.quantize(Operator::NuclearAttraction(molecule));
        ScalarSQOneElectronOperator<double> H = T + V;

        const auto g = g_spinor_basis.quantize(Operator::Coulomb());

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
    static enable_if_t<std::is_same<Z, double>::value, SQHamiltonian<double>> Random(const size_t K) {

        ScalarSQOneElectronOperator<double> H {SquareMatrix<double>::Random(K)};  // uniformly distributed between [-1,1]


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
     *  @param fcidump_filename     the name of the FCIDUMP file
     *
     *  @return the Hamiltonian corresponding to the contents of an FCIDUMP file
     *
     *  Note that this named constructor is only available for real matrix representations
     */
    template <typename Z = Scalar>
    static enable_if_t<std::is_same<Z, double>::value, SQHamiltonian<double>> ReadFCIDUMP(const std::string& fcidump_filename) {

        std::ifstream input_file_stream = validateAndOpen(fcidump_filename, "FCIDUMP");


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


        SquareMatrix<double> h_core = SquareMatrix<double>::Zero(K);
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
     *  Calculate the Edmiston-Ruedenberg localization index, which is the trace of the two-electron integrals over only the occupied orbitals.
     * 
     *  @param orbital_space                an orbital space which denotes the occupied-active-virtual separation
     *
     *  @return the Edmiston-Ruedenberg localization index
     *
     *  @note This method is only available for real matrix representations.
     */
    Scalar calculateEdmistonRuedenbergLocalizationIndex(const OrbitalSpace orbital_space) const {

        const auto& g_par = this->total_two_op.parameters();

        // TODO: when Eigen releases TensorTrace, use it here
        double localization_index = 0.0;
        for (const auto& i : orbital_space.indices(OccupationType::k_occupied)) {
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
     *  @param D            the 1-DM
     *  @param d            the 2-DM
     *
     *  @return the expectation value of this Hamiltonian
     */
    Scalar calculateExpectationValue(const OneDM<Scalar>& D, const TwoDM<Scalar>& d) const {

        return this->core().calculateExpectationValue(D)[0] + this->twoElectron().calculateExpectationValue(d)[0];  // SQHamiltonian contains ScalarSQOperators, so we access with [0]
    }


    /**
     *  @param D      the 1-DM (or the response 1-DM for made-variational wave function models)
     *  @param d      the 2-DM (or the response 2-DM for made-variational wave function models)
     *
     *  @return the (generalized) Fockian matrix
     */
    SquareMatrix<Scalar> calculateFockianMatrix(const OneDM<double>& D, const TwoDM<double>& d) const {

        return this->core().calculateFockianMatrix(D, d)[0] + this->twoElectron().calculateFockianMatrix(D, d)[0];  // SQHamiltonian has one- and two-electron contributions, so access with [0] accordingly
    }


    /**
     *  Calculate the (general) inactive Fockian operator.
     * 
     *  @param orbital_space                an orbital space which denotes the occupied-virtual separation
     * 
     *  @return the inactive Fockian operator
     */
    ScalarSQOneElectronOperator<Scalar> calculateInactiveFockian(const OrbitalSpace orbital_space) const {

        const auto& h_par = this->core().parameters();
        const auto& g_par = this->twoElectron().parameters();


        // A KISS implementation of the calculation of the (general) inactive Fockian matrix
        auto F_par = h_par;  // one-electron part

        // Two-electron part
        for (const auto& p : orbital_space.indices()) {
            for (const auto& q : orbital_space.indices()) {

                for (const auto& i : orbital_space.indices(OccupationType::k_occupied)) {
                    F_par(p, q) += g_par(p, q, i, i) - g_par(p, i, i, q);
                }
            }
        }  // F elements loop

        return ScalarSQOneElectronOperator<Scalar>(F_par);
    }


    /**
     *  Calculate the (restricted) inactive Fockian operator.
     * 
     *  @param orbital_space                an orbital space which denotes the occupied-virtual separation
     * 
     *  @return the inactive Fockian operator
     */
    ScalarSQOneElectronOperator<Scalar> calculateInactiveFockianRestricted(const OrbitalSpace orbital_space) const {

        const auto& h_par = this->core().parameters();
        const auto& g_par = this->twoElectron().parameters();


        // A KISS implementation of the calculation of the (restricted) inactive Fockian matrix
        auto F_par = h_par;  // one-electron part

        // Two-electron part
        for (const auto& p : orbital_space.indices()) {
            for (const auto& q : orbital_space.indices()) {

                for (const auto& i : orbital_space.indices(OccupationType::k_occupied)) {
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
    SquareRankFourTensor<Scalar> calculateSuperFockianMatrix(const OneDM<double>& D, const TwoDM<double>& d) const {

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
     *  @return the number of orbitals (spinors or spin-orbitals, depending on the context) this second-quantized Hamiltonian is related to
     */
    size_t numberOfOrbitals() const { return this->core().numberOfOrbitals(); }

    /**
     *  Using a random rotation matrix, transform the matrix representations of the Hamiltonian
     *
     *  Note that this method is only available for real matrix representations
     */
    template <typename Z = Scalar>
    enable_if_t<std::is_same<Z, double>::value> randomRotate() { this->rotate(TransformationMatrix<double>::RandomUnitary(this->numberOfOrbitals())); }

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


#include <numeric>
#include <vector>


/*
 *  MARK: SQHamiltonian
 */

/**
 *  A second-quantized electronic Hamiltonian. It consists of one-electron (core) and two-electron contributions.
 * 
 *  @tparam _ScalarSQOneElectronOperator        The type of second-quantized one-electron operator underlying this Hamiltonian.
 *  @tparam _ScalarSQTwoElectronOperator        The type of second-quantized two-electron operator underlying this Hamiltonian.
 */
template <typename _ScalarSQOneElectronOperator, typename _ScalarSQTwoElectronOperator>
class SQHamiltonian_Placeholder {
public:
    // The type of second-quantized one-electron operator underlying this Hamiltonian.
    using ScalarSQOneElectronOperator_Placeholder = _ScalarSQOneElectronOperator;  // TODO: remove 'placeholder'

    // The type of second-quantized two-electron operator underlying this Hamiltonian.
    using ScalarSQTwoElectronOperator_Placeholder = _ScalarSQTwoElectronOperator;  // TODO: remove 'placeholder'

    // Check if the spinor tags of the one- and two-electron operators match.
    assert(std::is_same<typename ScalarSQOneElectronOperator_Placeholder::SpinorTag, typename ScalarSQTwoElectronOperator_Placeholder::SpinorTag>::value, "The spinor tags of the one- and two-electron operators do not match.");

    // The spinor tag associated to this Hamiltonian.
    using SpinorTag = ScalarSQOneElectronOperator_Placeholder::SpinorTag;

    // The type of 'this'
    using Self = SQHamiltonian_Placeholder<ScalarSQOneElectronOperator_Placeholder, ScalarSQTwoElectronOperator_Placeholder>;


private:
    // The total one-electron interaction operator, i.e. the core Hamiltonian.
    ScalarSQOneElectronOperator_Placeholder h;

    // The total two-electron interaction operator.
    ScalarSQTwoElectronOperator_Placeholder g;

    // The contributions to the total one-electron interaction operator.
    std::vector<ScalarSQOneElectronOperator_Placeholder> h_contributions;

    // The contributions to the total two-electron interaction operator.
    std::vector<ScalarSQTwoElectronOperator_Placeholder> g_contributions;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  Initialize a `SQHamiltonian` for one- and two-electron contributions.
     * 
     *  @param h_contributions      The contributions to the total one-electron interaction operator.
     *  @param g_contributions      The contributions to the total two-electron interaction operator.
     */
    SQHamiltonian_Placeholder(const std::vector<ScalarSQOneElectronOperator_Placeholder>& h_contributions, const std::vector<ScalarSQTwoElectronOperator_Placeholder>& g_contributions) :
        h_contributions {h_contributions},
        g_contributions {g_contributions} {

        // Check if the dimensions of the one- and two-electron operators are compatible.
        const std::invalid_argument dimension_error {"SQHamiltonian::SQHamiltonian(const std::vector<ScalarSQOneElectronOperator<Scalar>& h_contributions, const std::vector<ScalarSQTwoElectronOperator<Scalar>& g_contributions: The dimensions of the contributing operators are incompatible"};

        const auto dim = h_contributions[0].numberOfOrbitals();
        for (const auto& h : this->h_contributions) {
            if (h.numberOfOrbitals() != dim) {
                throw dimension_error;
            }
        }

        for (const auto& g : this->g_contributions) {
            if (g.numberOfOrbitals() != dim) {
                throw dimension_error;
            }
        }


        // Calculate the total one- and two-electron operators.
        this->h = std::accumulate(h_contributions.begin(), h_contributions.end(), ScalarSQOneElectronOperator_Placeholder::Zero(dim));
        this->g = std::accumulate(g_contributions.begin(), g_contributions.end(), ScalarSQTwoElectronOperator_Placeholder::Zero(dim));
    }


    /**
     *  Create an `SQHamiltonian` from total one- and two-electron contributions.
     * 
     *  @param h            The total one-electron interaction operator, i.e. the core Hamiltonian.
     *  @param g            The total two-electron interaction operator.
     */
    SQHamiltonian_Placeholder(const ScalarSQOneElectronOperator_Placeholder& h, const ScalarSQTwoElectronOperator_Placeholder& g) :
        SQHamiltonian_Placeholder(std::vector<ScalarSQOneElectronOperator_Placeholder> {h},
                                  std::vector<ScalarSQTwoElectronOperator_Placeholder> {g}) {}


    /**
     *  MARK: Named constructors
     */

    /**
     *  Create an `SQHamiltonian` from a `HubbardHamiltonian`.
     * 
     *  @param hubbard_hamiltonian              The Hubbard model Hamiltonian.
     *
     *  @return An SQHamiltonian generated from the given Hubbard model Hamiltonian.
     *
     *  @note This named constructor is only available for real matrix representations.
     */
    template <typename Z1 = Scalar, typename Z2 = SpinorTag>
    static enable_if_t<std::is_same<Z1, double>::value && std::is_same<Z2, RestrictedSpinorTag>, SQHamiltonian_Placeholder<double>> FromHubbard(const GQCP::HubbardHamiltonian<double>& hubbard_hamiltonian) {

        const auto h = hubbard_hamiltonian.core();
        const auto g = hubbard_hamiltonian.twoElectron();

        return Self {h, g};
    }


    /**
     *  Construct the electronic molecular Hamiltonian in a spinor basis basis. The molecular Hamiltonian encompasses the following interactions:
     *      - one-electron contributions:
     *          - kinetic
     *          - nuclear attraction
     *      - two-electron contributions:
     *          - Coulomb repulsion
     * 
     *  @param spinor           The spinor basis in which the Hamiltonian should be expressed.
     *  @param molecule         The molecule that contains the nuclear framework upon which the nuclear attraction operator is based
     *
     *  @return An  second-quantized molecular Hamiltonian. 

     *
     *  @note This named constructor is only available for real matrix representations.
     */
    template <typename Scalar, typename SpinorBasis>
    static SQHamiltonian<typename SpinorBasis::ScalarSQOneElectronOperator_Placeholder, typename SpinorBasis::ScalarSQTwoElectronOperator_Placeholder> Molecular(const RSpinorBasis<Scalar, GTOShell>& spinor_basis, const Molecule& molecule) {

        // Calculate the integrals for the molecular Hamiltonian
        const auto T = spinor_basis.quantize(Operator::Kinetic());
        const auto V = spinor_basis.quantize(Operator::NuclearAttraction(molecule));
        const auto H = T + V;

        const auto g = spinor_basis.quantize(Operator::Coulomb());

        using ReturnType = SQHamiltonian<typename SpinorBasis::ScalarSQOneElectronOperator_Placeholder, typename SpinorBasis::ScalarSQTwoElectronOperator_Placeholder>;
        return ReturnType {H, g};
    }


    /**
     *  Create a random `SQHamiltonian` with values uniformly distributed between [-1,1].
     * 
     *  @param dim      The dimension of the matrix representations of the one- and two-electron integrals, i.e. the number of orbitals.
     *
     *  @return A random `SQHamiltonian`.
     *
     *  @note This named constructor is only available in the real case.
     */
    static enable_if_t<std::is_same<Z, double>::value, Self> Random(const size_t dim) {

        const auto h = ScalarSQOneElectronOperator_Placeholder::Random(dim);
        const auto g = ScalarSQTwoElectronOperator_Placeholder::Random(dim);

        return Self {h, g};
    }


    /**
     *  Parse an FCIDUMP file and read in its Hamiltonian.
     * 
     *  @param fcidump_filename         The name of the FCIDUMP file.
     *
     *  @return The Hamiltonian corresponding to the contents of an FCIDUMP file.
     *
     *  @note This named constructor is only available in the real case.
     */
    template <typename Z1 = Scalar, typename Z2 = SpinorTag>
    static enable_if_t<std::is_same<Z1, double>::value && std::is_same<Z2, RestrictedSpinOrbitalTag>::value, RSQHamiltonian<double>> FromFCIDUMP(const std::string& fcidump_filename) {

        std::ifstream input_file_stream = validateAndOpen(fcidump_filename, "FCIDUMP");


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


        SquareMatrix<double> h_core = SquareMatrix<double>::Zero(K);
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
     *  MARK: Calculations
     */

    /**
     *  Calculate the Edmiston-Ruedenberg localization index, which is the trace of the two-electron integrals over only the occupied orbitals.
     * 
     *  @param orbital_space                An orbital space which denotes the occupied-active-virtual separation.
     *
     *  @return The Edmiston-Ruedenberg localization index.
     *
     *  @note This method is only available for real matrix representations.
     */
    double calculateEdmistonRuedenbergLocalizationIndex(const OrbitalSpace orbital_space) const {

        const auto& g = this->g_total.parameters();

        // TODO: When Eigen releases TensorTrace, use it here.
        double localization_index = 0.0;
        for (const auto& i : orbital_space.indices(OccupationType::k_occupied)) {
            localization_index += g(i, i, i, i);
        }

        return localization_index;
    }


    /**
     *  Calculate the effective one-electron integrals. These are the core integrals, with the contributions form a Kronecker delta-term in the two-electron integrals.
     * 
     *  @return The effective one-electron integrals.
     */
    ScalarSQOneElectronOperator_Placeholder calculateEffectiveOneElectronIntegrals() const {

        return this->core() + this->twoElectron().effectiveOneElectronPartition();
    }
};


/*
 *  MARK: Convenience aliases
 */

// An `SQHamiltonian` related to restricted spin-orbitals. See `RestrictedSpinOrbitalTag`.
template <typename Scalar>
using RSQHamiltonian = SQHamiltonian_Placeholder<ScalarRSQOneElectronOperator<Scalar>, ScalarSQTwoElectronOperator<Scalar>>;

// An `SQHamiltonian` related to unrestricted spin-orbitals. See `UnrestrictedSpinOrbitalTag`.
// template <typename Scalar>
// using USQHamiltonian = SQHamiltonian_Placeholder<ScalarUSQOneElectronOperator<Scalar>, ScalarUSQTwoElectronOperator<Scalar>>;

// An `SQHamiltonian` related to general spinors. See `GeneralSpinorTag`.
template <typename Scalar>
using GSQHamiltonian = GQHamiltonian_Placeholder<ScalarGSQOneElectronOperator<Scalar>, ScalarSQTwoElectronOperator<Scalar>>;


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
