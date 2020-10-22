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


#include "Basis/SpinorBasis/GSpinorBasis.hpp"
#include "Basis/SpinorBasis/OrbitalSpace.hpp"
#include "Basis/SpinorBasis/RSpinorBasis.hpp"
#include "Basis/Transformations/BasisTransformable.hpp"
#include "Basis/Transformations/JacobiRotatable.hpp"
#include "Operator/SecondQuantized/GSQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/GSQTwoElectronOperator.hpp"
#include "Operator/SecondQuantized/ModelHamiltonian/HubbardHamiltonian.hpp"
#include "Operator/SecondQuantized/RSQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/RSQTwoElectronOperator.hpp"
#include "QuantumChemical/spinor_tags.hpp"
#include "Utilities/type_traits.hpp"

#include <numeric>
#include <vector>


namespace GQCP {


/*
 *  MARK: SQHamiltonian implementation
 */

/**
 *  A second-quantized electronic Hamiltonian. It consists of one-electron (core) and two-electron contributions.
 * 
 *  @tparam _ScalarSQOneElectronOperator        The type of second-quantized one-electron operator underlying this Hamiltonian.
 *  @tparam _ScalarSQTwoElectronOperator        The type of second-quantized two-electron operator underlying this Hamiltonian.
 */
template <typename ScalarSQOneElectronOperator_Placeholder, typename ScalarSQTwoElectronOperator_Placeholder>
class SQHamiltonian:
    public BasisTransformable<SQHamiltonian<ScalarSQOneElectronOperator_Placeholder, ScalarSQTwoElectronOperator_Placeholder>>,
    public JacobiRotatable<SQHamiltonian<ScalarSQOneElectronOperator_Placeholder, ScalarSQTwoElectronOperator_Placeholder>> {
public:
    // // The type of second-quantized one-electron operator underlying this Hamiltonian.
    // using ScalarSQOneElectronOperator_Placeholder = ScalarSQOneElectronOperator_Placeholder;  // TODO: remove 'placeholder'

    // // The type of second-quantized two-electron operator underlying this Hamiltonian.
    // using ScalarSQTwoElectronOperator_Placeholder = ScalarSQTwoElectronOperator_Placeholder;  // TODO: remove 'placeholder'

    // Check if the spinor tags of the one- and two-electron operators match.
    static_assert(std::is_same<typename ScalarSQOneElectronOperator_Placeholder::SpinorTag, typename ScalarSQTwoElectronOperator_Placeholder::SpinorTag>::value, "The spinor tags of the one- and two-electron operators do not match.");

    // The spinor tag associated to this Hamiltonian.
    using SpinorTag = typename ScalarSQOneElectronOperator_Placeholder::SpinorTag;

    // Check if the scalar type of the parameters/matrix elements of the one-electron and two-electron operators are equal.
    static_assert(std::is_same<typename ScalarSQOneElectronOperator_Placeholder::Scalar, typename ScalarSQTwoElectronOperator_Placeholder::Scalar>::value, "The scalar type of the one- and two-electron parameters/matrix elements do not match.");

    // The scalar type used for a single parameter/matrix element: real or complex.
    using Scalar = typename ScalarSQOneElectronOperator_Placeholder::Scalar;

    // The type of 'this'
    using Self = SQHamiltonian<ScalarSQOneElectronOperator_Placeholder, ScalarSQTwoElectronOperator_Placeholder>;

    // The type of the one-particle density matrix that is naturally associated to the second-quantized Hamiltonian.
    using OneDM = typename OperatorTraits<ScalarSQOneElectronOperator_Placeholder>::OneDM;

    // The type of the two-particle density matrix that is naturally associated to the second-quantized Hamiltonian.
    using TwoDM = typename OperatorTraits<ScalarSQOneElectronOperator_Placeholder>::TwoDM;

    // The type of transformation matrix that is naturally associated to the Hamiltonian.
    using TM = typename OperatorTraits<ScalarSQOneElectronOperator_Placeholder>::TM;


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
    SQHamiltonian(const std::vector<ScalarSQOneElectronOperator_Placeholder>& h_contributions, const std::vector<ScalarSQTwoElectronOperator_Placeholder>& g_contributions) :
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
    SQHamiltonian(const ScalarSQOneElectronOperator_Placeholder& h, const ScalarSQTwoElectronOperator_Placeholder& g) :
        SQHamiltonian(std::vector<ScalarSQOneElectronOperator_Placeholder> {h},
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
    static enable_if_t<std::is_same<Z1, double>::value && std::is_same<Z2, RestrictedSpinorTag>::value, SQHamiltonian<ScalarSQOneElectronOperator_Placeholder, ScalarSQTwoElectronOperator_Placeholder>> FromHubbard(const GQCP::HubbardHamiltonian<double>& hubbard_hamiltonian) {

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
     *  @return A second-quantized molecular Hamiltonian.
     */
    // FIXME: This API should be moved to SpinorBasis.
    template <typename Z = SpinorTag>
    static enable_if_t<std::is_same<Z, RestrictedSpinOrbitalTag>::value, SQHamiltonian<ScalarRSQOneElectronOperator<double>, ScalarRSQTwoElectronOperator<double>>> Molecular(const RSpinorBasis<double, GTOShell>& spinor_basis, const Molecule& molecule) {

        // Calculate the integrals for the molecular Hamiltonian
        const auto T = spinor_basis.quantize(Operator::Kinetic());
        const auto V = spinor_basis.quantize(Operator::NuclearAttraction(molecule));
        const auto H = T + V;

        const auto g = spinor_basis.quantize(Operator::Coulomb());

        return SQHamiltonian<ScalarRSQOneElectronOperator<double>, ScalarRSQTwoElectronOperator<double>> {H, g};
    }
    template <typename Z = SpinorTag>
    static enable_if_t<std::is_same<Z, GeneralSpinorTag>::value, SQHamiltonian<ScalarGSQOneElectronOperator<double>, ScalarGSQTwoElectronOperator<double>>> Molecular(const GSpinorBasis<double, GTOShell>& spinor_basis, const Molecule& molecule) {

        // Calculate the integrals for the molecular Hamiltonian
        const auto T = spinor_basis.quantize(Operator::Kinetic());
        const auto V = spinor_basis.quantize(Operator::NuclearAttraction(molecule));
        const auto H = T + V;

        const auto g = spinor_basis.quantize(Operator::Coulomb());

        return SQHamiltonian<ScalarGSQOneElectronOperator<double>, ScalarGSQTwoElectronOperator<double>> {H, g};
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
    template <typename Z = Scalar>
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
    static enable_if_t<std::is_same<Z1, double>::value && std::is_same<Z2, RestrictedSpinOrbitalTag>::value, SQHamiltonian<ScalarSQOneElectronOperator_Placeholder, ScalarSQTwoElectronOperator_Placeholder>> FromFCIDUMP(const std::string& fcidump_filename) {

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
            throw std::invalid_argument("SQHamiltonian::FromFCIDUMP(std::string): The .FCIDUMP-file is invalid: could not read a number of orbitals.");
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


        return SQHamiltonian(ScalarSQOneElectronOperator_Placeholder(h_core), ScalarSQTwoElectronOperator_Placeholder(g));
    }


    /*
     *  MARK: Parameter access
     */

    /**
     *  @return A read-only reference to the total one-electron interaction operator, i.e. the core Hamiltonian.
     */
    const ScalarSQOneElectronOperator_Placeholder& core() const { return this->h; }

    /**
     *  @return A writable reference to the total one-electron interaction operator, i.e. the core Hamiltonian.
     */
    ScalarSQOneElectronOperator_Placeholder& core() { return this->h; }

    /**
     *  @return A read-only reference to the contributions to the total one-electron interaction operator.
     */
    const std::vector<ScalarSQOneElectronOperator_Placeholder>& coreContributions() const { return this->h_contributions; }

    /**
     *  @return A writable reference to the contributions to the total one-electron interaction operator.
     */
    std::vector<ScalarSQOneElectronOperator_Placeholder>& coreContributions() { return this->h_contributions; }

    /**
     *  @return A read-only reference to the total two-electron interaction operator.
     */
    const ScalarSQTwoElectronOperator_Placeholder& twoElectron() const { return this->g; }

    /**
     *  @return A writable reference to the total two-electron interaction operator.
     */
    ScalarSQTwoElectronOperator_Placeholder& twoElectron() { return this->g; }

    /**
     *  @return A read-only reference to the contributions to the total two-electron interaction operator.
     */
    const std::vector<ScalarSQTwoElectronOperator_Placeholder>& twoElectronContributions() const { return this->g_contributions; }

    /**
     *  @return A writable reference to the contributions to the total two-electron interaction operator.
     */
    std::vector<ScalarSQTwoElectronOperator_Placeholder>& twoElectronContributions() { return this->g_contributions; }


    /*
     *  MARK: General information
     */

    /**
     *  @return The number of orbitals (spinors or spin-orbitals, depending on the context) this second-quantized Hamiltonian is related to
     */
    size_t numberOfOrbitals() const { return this->core().numberOfOrbitals(); }


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

        const auto& g_total_par = this->twoElectron().parameters();

        // TODO: When Eigen releases TensorTrace, use it here.
        double localization_index = 0.0;
        for (const auto& i : orbital_space.indices(OccupationType::k_occupied)) {
            localization_index += g_total_par(i, i, i, i);
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


    /**
     *  Calculate the electronic energy, i.e. the expectation value of this Hamiltonian.
     * 
     *  @param D            The 1-DM.
     *  @param d            The 2-DM.
     *
     *  @return The expectation value of this Hamiltonian.
     */
    Scalar calculateExpectationValue(const OneDM& D, const TwoDM& d) const {

        // An SQHamiltonian contains ScalarSQOperators, so we access their expectation values with ().
        return this->core().calculateExpectationValue(D)() + this->twoElectron().calculateExpectationValue(d)();
    }


    /**
     *  Calculate the Fockian matrix of this Hamiltonian.
     * 
     *  @param D      The 1-DM (or the response 1-DM for made-variational wave function models).
     *  @param d      The 2-DM (or the response 2-DM for made-variational wave function models).
     *
     *  @return The Fockian matrix.
     */
    SquareMatrix<Scalar> calculateFockianMatrix(const OneDM& D, const TwoDM& d) const {

        // An SQHamiltonian contains ScalarSQOperators, so we access their Fockian matrices with (0).
        return this->core().calculateFockianMatrix(D, d)() + this->twoElectron().calculateFockianMatrix(D, d)();
    }


    /**
     *  Calculate the super-Fockian matrix of this Hamiltonian.
     * 
     *  @param D      The 1-DM (or the response 1-DM for made-variational wave function models).
     *  @param d      The 2-DM (or the response 2-DM for made-variational wave function models).
     *
     *  @return The super-Fockian matrix.
     */
    SquareRankFourTensor<Scalar> calculateSuperFockianMatrix(const OneDM& D, const TwoDM& d) const {

        // An SQHamiltonian contains ScalarSQOperators, so we access their Fockian matrices with (0).
        return this->core().calculateSuperFockianMatrix(D, d)().Eigen() + this->twoElectron().calculateSuperFockianMatrix(D, d)().Eigen();  // We have to call .Eigen() because operator+ isn't enabled on SquareRankFourTensor.
    }


    /**
     *  Calculate the inactive Fockian operator.
     * 
     *  @param orbital_space                An orbital space which denotes the occupied-virtual separation.
     * 
     *  @return The inactive Fockian operator.
     */
    template <typename Z = SpinorTag>
    enable_if_t<std::is_same<Z, GeneralSpinorTag>::value, ScalarSQOneElectronOperator_Placeholder> calculateInactiveFockian(const OrbitalSpace orbital_space) const {

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

        return ScalarSQOneElectronOperator_Placeholder {F_par};
    }


    /**
     *  Calculate the inactive Fockian operator.
     * 
     *  @param orbital_space                An orbital space which denotes the occupied-virtual separation.
     * 
     *  @return The inactive Fockian operator.
     */
    template <typename Z = SpinorTag>
    enable_if_t<std::is_same<Z, RestrictedSpinOrbitalTag>::value, ScalarSQOneElectronOperator_Placeholder> calculateInactiveFockian(const OrbitalSpace orbital_space) const {
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

        return ScalarSQOneElectronOperator_Placeholder {F_par};
    }


    /*
     *  MARK: Conforming to BasisTransformable
     */

    /**
     *  Apply the basis transformation and return the resulting one-electron integrals.
     * 
     *  @param T            The type that encapsulates the basis transformation coefficients.
     * 
     *  @return The basis-transformed one-electron integrals.
     */
    Self transformed(const TM& T) const override {

        auto result = *this;

        // Transform the one and two-electron contributions.
        for (auto& h : result.coreContributions()) {
            h.transform(T);
        }

        for (auto& g : result.twoElectronContributions()) {
            g.transform(T);
        }

        // Transform the total one- and two-electron interactions.
        result.core().transform(T);
        result.twoElectron().transform(T);

        return result;
    }


    // Allow the `rotate` method from `BasisTransformable`, since there's also a `rotate` from `JacobiRotatable`.
    using BasisTransformable<Self>::rotate;

    // Allow the `rotated` method from `BasisTransformable`, since there's also a `rotated` from `JacobiRotatable`.
    using BasisTransformable<Self>::rotated;


    /*
     *  MARK: Conforming to JacobiRotatable
     */

    /**
     *  Apply the Jacobi rotation and return the result.
     * 
     *  @param jacobi_parameters        The Jacobi rotation parameters.
     * 
     *  @return The jacobi-transformed object.
     */
    Self rotated(const JacobiRotationParameters& jacobi_parameters) const override {

        auto result = *this;

        // Transform the one and two-electron contributions.
        for (auto& h : result.coreContributions()) {
            h.rotate(jacobi_parameters);
        }

        for (auto& g : result.twoElectronContributions()) {
            g.rotate(jacobi_parameters);
        }

        // Transform the total one- and two-electron interactions.
        result.core().rotate(jacobi_parameters);
        result.twoElectron().rotate(jacobi_parameters);

        return result;
    }

    // Allow the `rotate` method from `JacobiRotatable`, since there's also a `rotate` from `BasisTransformable`.
    using JacobiRotatable<Self>::rotate;
};


/*
 *  MARK: Operator traits
 */

/**
 *  A type that provides compile-time information on operators that is otherwise not accessible through a public class alias.
 */
template <typename ScalarSQOneElectronOperator_Placeholder, typename ScalarSQTwoElectronOperator_Placeholder>
struct OperatorTraits<SQHamiltonian<ScalarSQOneElectronOperator_Placeholder, ScalarSQTwoElectronOperator_Placeholder>> {

    // The type of transformation matrix that is naturally associated to the second-quantized Hamiltonian.
    using TM = typename OperatorTraits<ScalarSQOneElectronOperator_Placeholder>::TM;

    // The type of the one-particle density matrix that is naturally associated to the second-quantized Hamiltonian.
    using OneDM = typename OperatorTraits<ScalarSQOneElectronOperator_Placeholder>::OneDM;

    // The type of the two-particle density matrix that is naturally associated to the second-quantized Hamiltonian.
    using TwoDM = typename OperatorTraits<ScalarSQOneElectronOperator_Placeholder>::TwoDM;
};


/*
 *  MARK: BasisTransformableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `BasisTransformable`.
 */
template <typename ScalarSQOneElectronOperator_Placeholder, typename ScalarSQTwoElectronOperator_Placeholder>
struct BasisTransformableTraits<SQHamiltonian<ScalarSQOneElectronOperator_Placeholder, ScalarSQTwoElectronOperator_Placeholder>> {

    // The type of the transformation matrix for which the basis transformation should be defined. // TODO: Rename "TM" to "TransformationMatrix"
    using TM = typename OperatorTraits<ScalarSQOneElectronOperator_Placeholder>::TM;

    // using TM = typename OperatorTraits<SQHamiltonian<ScalarSQOneElectronOperator_Placeholder, ScalarSQTwoElectronOperator_Placeholder>>::TM;
};


/*
 *  MARK: Convenience aliases
 */


// An `SQHamiltonian` related to restricted spin-orbitals. See `RestrictedSpinOrbitalTag`.
template <typename Scalar>
using RSQHamiltonian = SQHamiltonian<ScalarRSQOneElectronOperator<Scalar>, ScalarRSQTwoElectronOperator<Scalar>>;


// An `SQHamiltonian` related to unrestricted spin-orbitals. See `UnrestrictedSpinOrbitalTag`.
// template <typename Scalar>
// using USQHamiltonian = SQHamiltonian<ScalarUSQOneElectronOperator<Scalar>, ScalarUSQTwoElectronOperator<Scalar>>;


// An `SQHamiltonian` related to general spinors. See `GeneralSpinorTag`.
template <typename Scalar>
using GSQHamiltonian = SQHamiltonian<ScalarGSQOneElectronOperator<Scalar>, ScalarGSQTwoElectronOperator<Scalar>>;


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
RSQHamiltonian<Scalar> operator+(const RSQHamiltonian<Scalar>& sq_hamiltonian, const ScalarRSQOneElectronOperator<Scalar>& sq_one_op) {

    // Make a copy of the one-electron part in order to create a new Hamiltonian
    auto sq_one_ops = sq_hamiltonian.coreContributions();

    // 'Add' the one-electron operator
    sq_one_ops.push_back(sq_one_op);

    return RSQHamiltonian<Scalar>(sq_one_ops, sq_hamiltonian.twoElectronContributions());
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
RSQHamiltonian<Scalar> operator-(const RSQHamiltonian<Scalar>& sq_hamiltonian, const ScalarRSQOneElectronOperator<Scalar>& sq_one_op) {

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
RSQHamiltonian<Scalar> operator+(const RSQHamiltonian<Scalar>& sq_hamiltonian, const ScalarRSQTwoElectronOperator<Scalar>& sq_two_op) {

    // Make a copy of the two-electron part in order to create a new Hamiltonian
    auto sq_two_ops = sq_hamiltonian.twoElectronContributions();

    // 'Add' the two-electron operator
    sq_two_ops.push_back(sq_two_ops);

    return RSQHamiltonian<Scalar>(sq_hamiltonian.coreContributions(), sq_two_ops);
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
RSQHamiltonian<Scalar> operator-(const RSQHamiltonian<Scalar>& sq_hamiltonian, const ScalarRSQTwoElectronOperator<Scalar>& sq_two_op) {

    return sq_hamiltonian + (-sq_two_op);
}


}  // namespace GQCP
