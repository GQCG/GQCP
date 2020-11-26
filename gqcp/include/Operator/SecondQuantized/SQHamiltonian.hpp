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
#include "Basis/SpinorBasis/RSpinOrbitalBasis.hpp"
#include "Basis/SpinorBasis/USpinOrbitalBasis.hpp"
#include "Basis/Transformations/BasisTransformable.hpp"
#include "Basis/Transformations/JacobiRotatable.hpp"
#include "Operator/SecondQuantized/GSQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/GSQTwoElectronOperator.hpp"
#include "Operator/SecondQuantized/ModelHamiltonian/HubbardHamiltonian.hpp"
#include "Operator/SecondQuantized/RSQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/RSQTwoElectronOperator.hpp"
#include "Operator/SecondQuantized/USQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/USQTwoElectronOperator.hpp"
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
template <typename _ScalarSQOneElectronOperator, typename _ScalarSQTwoElectronOperator>
class SQHamiltonian:
    public BasisTransformable<SQHamiltonian<_ScalarSQOneElectronOperator, _ScalarSQTwoElectronOperator>>,
    public JacobiRotatable<SQHamiltonian<_ScalarSQOneElectronOperator, _ScalarSQTwoElectronOperator>> {
public:
    // The type of second-quantized one-electron operator underlying this Hamiltonian.
    using ScalarSQOneElectronOperator = _ScalarSQOneElectronOperator;

    // The type of second-quantized two-electron operator underlying this Hamiltonian.
    using ScalarSQTwoElectronOperator = _ScalarSQTwoElectronOperator;

    // Check if the spinor tags of the one- and two-electron operators match.
    static_assert(std::is_same<typename ScalarSQOneElectronOperator::SpinorTag, typename ScalarSQTwoElectronOperator::SpinorTag>::value, "The spinor tags of the one- and two-electron operators do not match.");

    // The spinor tag associated to this Hamiltonian.
    using SpinorTag = typename ScalarSQOneElectronOperator::SpinorTag;

    // Check if the scalar type of the parameters/matrix elements of the one-electron and two-electron operators are equal.
    static_assert(std::is_same<typename ScalarSQOneElectronOperator::Scalar, typename ScalarSQTwoElectronOperator::Scalar>::value, "The scalar type of the one- and two-electron parameters/matrix elements do not match.");

    // The scalar type used for a single parameter/matrix element: real or complex.
    using Scalar = typename ScalarSQOneElectronOperator::Scalar;

    // The type of 'this'
    using Self = SQHamiltonian<ScalarSQOneElectronOperator, ScalarSQTwoElectronOperator>;

    // The type of the one-particle density matrix that is naturally associated to the second-quantized Hamiltonian.
    using OneDM = typename OperatorTraits<ScalarSQOneElectronOperator>::OneDM;

    // The type of the two-particle density matrix that is naturally associated to the second-quantized Hamiltonian.
    using TwoDM = typename OperatorTraits<ScalarSQOneElectronOperator>::TwoDM;

    // The type of transformation that is naturally associated to the Hamiltonian.
    using Transformation = typename OperatorTraits<ScalarSQOneElectronOperator>::Transformation;

    // The type of Jacobi rotation for which the Jacobi rotation should be defined.
    using JacobiRotationType = typename JacobiRotatableTraits<ScalarSQOneElectronOperator>::JacobiRotationType;


private:
    // The total one-electron interaction operator, i.e. the core Hamiltonian.
    ScalarSQOneElectronOperator h;

    // The total two-electron interaction operator.
    ScalarSQTwoElectronOperator g;

    // The contributions to the total one-electron interaction operator.
    std::vector<ScalarSQOneElectronOperator> h_contributions;

    // The contributions to the total two-electron interaction operator.
    std::vector<ScalarSQTwoElectronOperator> g_contributions;


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
    SQHamiltonian(const std::vector<ScalarSQOneElectronOperator>& h_contributions, const std::vector<ScalarSQTwoElectronOperator>& g_contributions) :
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
        this->h = std::accumulate(h_contributions.begin(), h_contributions.end(), ScalarSQOneElectronOperator::Zero(dim));
        this->g = std::accumulate(g_contributions.begin(), g_contributions.end(), ScalarSQTwoElectronOperator::Zero(dim));
    }


    /**
     *  Create an `SQHamiltonian` from total one- and two-electron contributions.
     * 
     *  @param h            The total one-electron interaction operator, i.e. the core Hamiltonian.
     *  @param g            The total two-electron interaction operator.
     */
    SQHamiltonian(const ScalarSQOneElectronOperator& h, const ScalarSQTwoElectronOperator& g) :
        SQHamiltonian(std::vector<ScalarSQOneElectronOperator> {h},
                      std::vector<ScalarSQTwoElectronOperator> {g}) {}


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
    static enable_if_t<std::is_same<Z1, double>::value && std::is_same<Z2, RestrictedSpinorTag>::value, SQHamiltonian<ScalarSQOneElectronOperator, ScalarSQTwoElectronOperator>> FromHubbard(const GQCP::HubbardHamiltonian<double>& hubbard_hamiltonian) {

        const auto h = hubbard_hamiltonian.core();
        const auto g = hubbard_hamiltonian.twoElectron();

        return Self {h, g};
    }


    // FIXME: These APIs should be moved to SpinorBasis.

    /**
     *  Quantize the molecular Hamiltonian in a restricted spin-orbital basis.
     * 
     *  @param spinor_basis         The spinor basis in which the Hamiltonian should be expressed.
     *  @param molecule             The molecule that contains the nuclear framework upon which the nuclear attraction operator is based.
     *
     *  @return A second-quantized molecular Hamiltonian.
     */
    template <typename Z = SpinorTag>
    static enable_if_t<std::is_same<Z, RestrictedSpinOrbitalTag>::value, Self> Molecular(const RSpinOrbitalBasis<double, GTOShell>& spinor_basis, const Molecule& molecule) {

        // Calculate the integrals for the molecular Hamiltonian
        const auto T = spinor_basis.quantize(Operator::Kinetic());
        const auto V = spinor_basis.quantize(Operator::NuclearAttraction(molecule));
        const auto H = T + V;

        const auto g = spinor_basis.quantize(Operator::Coulomb());

        return Self {H, g};
    }

    /**
     *  Quantize the molecular Hamiltonian in an unrestricted spin-orbital basis.
     * 
     *  @param spinor_basis         The spinor basis in which the Hamiltonian should be expressed.
     *  @param molecule             The molecule that contains the nuclear framework upon which the nuclear attraction operator is based.
     *
     *  @return A second-quantized molecular Hamiltonian.
     */
    template <typename Z = SpinorTag>
    static enable_if_t<std::is_same<Z, UnrestrictedSpinOrbitalTag>::value, Self> Molecular(const USpinOrbitalBasis<double, GTOShell>& spinor_basis, const Molecule& molecule) {

        // Calculate the integrals for the molecular Hamiltonian
        const auto T = spinor_basis.quantize(Operator::Kinetic());
        const auto V = spinor_basis.quantize(Operator::NuclearAttraction(molecule));
        const auto H = T + V;

        const auto g = spinor_basis.quantize(Operator::Coulomb());

        return Self {H, g};
    }

    /**
     *  Quantize the molecular Hamiltonian in a general spinor basis.
     * 
     *  @param spinor_basis         The spinor basis in which the Hamiltonian should be expressed.
     *  @param molecule             The molecule that contains the nuclear framework upon which the nuclear attraction operator is based.
     *
     *  @return A second-quantized molecular Hamiltonian.
     */
    template <typename Z = SpinorTag>
    static enable_if_t<std::is_same<Z, GeneralSpinorTag>::value, Self> Molecular(const GSpinorBasis<double, GTOShell>& spinor_basis, const Molecule& molecule) {

        // Calculate the integrals for the molecular Hamiltonian
        const auto T = spinor_basis.quantize(Operator::Kinetic());
        const auto V = spinor_basis.quantize(Operator::NuclearAttraction(molecule));
        const auto H = T + V;

        const auto g = spinor_basis.quantize(Operator::Coulomb());

        return Self {H, g};
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

        const auto h = ScalarSQOneElectronOperator::Random(dim);
        const auto g = ScalarSQTwoElectronOperator::Random(dim);

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
    static enable_if_t<std::is_same<Z1, double>::value && std::is_same<Z2, RestrictedSpinOrbitalTag>::value, SQHamiltonian<ScalarSQOneElectronOperator, ScalarSQTwoElectronOperator>> FromFCIDUMP(const std::string& fcidump_filename) {

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
        SquareRankFourTensor<double> g {K};
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


        return SQHamiltonian(ScalarSQOneElectronOperator(h_core), ScalarSQTwoElectronOperator(g));
    }


    /*
     *  MARK: Access
     */

    /**
     *  @return A read-only reference to the total one-electron interaction operator, i.e. the core Hamiltonian.
     */
    const ScalarSQOneElectronOperator& core() const { return this->h; }

    /**
     *  @return A writable reference to the total one-electron interaction operator, i.e. the core Hamiltonian.
     */
    ScalarSQOneElectronOperator& core() { return this->h; }

    /**
     *  @return A read-only reference to the contributions to the total one-electron interaction operator.
     */
    const std::vector<ScalarSQOneElectronOperator>& coreContributions() const { return this->h_contributions; }

    /**
     *  @return A writable reference to the contributions to the total one-electron interaction operator.
     */
    std::vector<ScalarSQOneElectronOperator>& coreContributions() { return this->h_contributions; }

    /**
     *  @return A read-only reference to the total two-electron interaction operator.
     */
    const ScalarSQTwoElectronOperator& twoElectron() const { return this->g; }

    /**
     *  @return A writable reference to the total two-electron interaction operator.
     */
    ScalarSQTwoElectronOperator& twoElectron() { return this->g; }

    /**
     *  @return A read-only reference to the contributions to the total two-electron interaction operator.
     */
    const std::vector<ScalarSQTwoElectronOperator>& twoElectronContributions() const { return this->g_contributions; }

    /**
     *  @return A writable reference to the contributions to the total two-electron interaction operator.
     */
    std::vector<ScalarSQTwoElectronOperator>& twoElectronContributions() { return this->g_contributions; }


    /*
     *  MARK: General information
     */

    /**
     *  @return The number of orbitals (spinors or spin-orbitals, depending on the context) this second-quantized Hamiltonian is related to.
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
     *  @note This method is not enabled for unrestricted Hamiltonians.
     */
    template <typename Z = SpinorTag>
    enable_if_t<std::is_same<Z, RestrictedSpinOrbitalTag>::value || std::is_same<Z, GeneralSpinorTag>::value, double> calculateEdmistonRuedenbergLocalizationIndex(const OrbitalSpace orbital_space) const {

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
     * 
     *  @note This method is not enabled for unrestricted Hamiltonians.
     */
    template <typename Z = SpinorTag>
    enable_if_t<std::is_same<Z, RestrictedSpinOrbitalTag>::value || std::is_same<Z, GeneralSpinorTag>::value, ScalarSQOneElectronOperator> calculateEffectiveOneElectronIntegrals() const {

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
     * 
     *  @note This method is not enabled for unrestricted Hamiltonians.
     */
    template <typename Z = SpinorTag>
    enable_if_t<std::is_same<Z, RestrictedSpinOrbitalTag>::value || std::is_same<Z, GeneralSpinorTag>::value, SquareMatrix<Scalar>> calculateFockianMatrix(const OneDM& D, const TwoDM& d) const {

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
    template <typename Z = SpinorTag>
    enable_if_t<std::is_same<Z, RestrictedSpinOrbitalTag>::value || std::is_same<Z, GeneralSpinorTag>::value, SquareRankFourTensor<Scalar>> calculateSuperFockianMatrix(const OneDM& D, const TwoDM& d) const {

        // An SQHamiltonian contains ScalarSQOperators, so we access their Fockian matrices with (0).
        return this->core().calculateSuperFockianMatrix(D, d)().Eigen() + this->twoElectron().calculateSuperFockianMatrix(D, d)().Eigen();  // We have to call .Eigen() because operator+ isn't fully enabled on SquareRankFourTensor.
    }


    /**
     *  Calculate the inactive Fockian operator for a general Hamiltonian.
     * 
     *  @param orbital_space                An orbital space which denotes the occupied-virtual separation.
     * 
     *  @return The inactive Fockian operator.
     */
    template <typename Z = SpinorTag>
    enable_if_t<std::is_same<Z, GeneralSpinorTag>::value, ScalarSQOneElectronOperator> calculateInactiveFockian(const OrbitalSpace orbital_space) const {

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

        return ScalarSQOneElectronOperator {F_par};
    }


    /**
     *  Calculate the inactive Fockian operator for a restricted Hamiltonian.
     * 
     *  @param orbital_space                An orbital space which denotes the occupied-virtual separation.
     * 
     *  @return The inactive Fockian operator.
     */
    template <typename Z = SpinorTag>
    enable_if_t<std::is_same<Z, RestrictedSpinOrbitalTag>::value, ScalarSQOneElectronOperator> calculateInactiveFockian(const OrbitalSpace orbital_space) const {
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

        return ScalarSQOneElectronOperator {F_par};
    }


    /*
     *  MARK: Conforming to BasisTransformable
     */

    /**
     *  Apply the basis transformation and return the resulting Hamiltonian.
     * 
     *  @param T            The basis transformation.
     * 
     *  @return The basis-transformed Hamiltonian.
     */
    Self transformed(const Transformation& T) const override {

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
     *  @param jacobi_rotation          The Jacobi rotation.
     * 
     *  @return The Jacobi-rotated object.
     */
    Self rotated(const JacobiRotationType& jacobi_rotation) const override {

        auto result = *this;

        // Transform the one and two-electron contributions.
        for (auto& h : result.coreContributions()) {
            h.rotate(jacobi_rotation);
        }

        for (auto& g : result.twoElectronContributions()) {
            g.rotate(jacobi_rotation);
        }

        // Transform the total one- and two-electron interactions.
        result.core().rotate(jacobi_rotation);
        result.twoElectron().rotate(jacobi_rotation);

        return result;
    }

    // Allow the `rotate` method from `JacobiRotatable`, since there's also a `rotate` from `BasisTransformable`.
    using JacobiRotatable<Self>::rotate;


    /*
     *  MARK: Operations related to one-electron operators
     */

    /**
     *  Addition-assignment with a (scalar) one-electron operator.
     */
    Self& operator+=(const ScalarSQOneElectronOperator& sq_one_op) {

        // Update the one-electron contributions and the total one-electron operator.
        this->h_contributions.push_back(sq_one_op);
        this->h += sq_one_op;

        return *this;
    }


    /**
     *  Addition, canonically implemented using addition-assignment.
     */
    friend Self operator+(Self lhs, const ScalarSQOneElectronOperator& rhs) {
        lhs += rhs;
        return lhs;
    }


    /**
     *  A commutative version of the previous addition.
     */
    friend Self operator+(const ScalarSQOneElectronOperator& lhs, Self rhs) {
        return rhs += lhs;
    }


    /**
     *  Subtraction-assignment with a (scalar) one-electron operator.
     */
    Self& operator-=(const ScalarSQOneElectronOperator& sq_one_op) {

        *this += -sq_one_op;
        return *this;
    }


    /**
     *  Subtraction, canonically implemented using addition-assignment.
     */
    friend Self operator-(Self lhs, const ScalarSQOneElectronOperator& rhs) {
        lhs -= rhs;
        return lhs;
    }


    /*
     *  MARK: Operations related to two-electron operators
     */

    /**
     *  Addition-assignment with a (scalar) two-electron operator.
     */
    Self& operator+=(const ScalarSQTwoElectronOperator& sq_two_op) {

        // Update the one-electron contributions and the total one-electron operator.
        this->g_contributions.push_back(sq_two_op);
        this->g += sq_two_op;

        return *this;
    }


    /**
     *  Addition, canonically implemented using addition-assignment.
     */
    friend Self operator+(Self lhs, const ScalarSQTwoElectronOperator& rhs) {
        lhs += rhs;
        return lhs;
    }


    /**
     *  A commutative version of the previous addition.
     */
    friend Self operator+(const ScalarSQTwoElectronOperator& lhs, Self rhs) {
        return rhs += lhs;
    }


    /**
     *  Subtraction-assignment with a (scalar) one-electron operator.
     */
    Self& operator-=(const ScalarSQTwoElectronOperator& sq_two_op) {

        *this += -sq_two_op;
        return *this;
    }


    /**
     *  Subtraction, canonically implemented using addition-assignment.
     */
    friend Self operator-(Self lhs, const ScalarSQTwoElectronOperator& rhs) {
        lhs -= rhs;
        return lhs;
    }
};


/*
 *  MARK: Operator traits
 */

/**
 *  A type that provides compile-time information on operators that is otherwise not accessible through a public class alias.
 */
template <typename ScalarSQOneElectronOperator, typename ScalarSQTwoElectronOperator>
struct OperatorTraits<SQHamiltonian<ScalarSQOneElectronOperator, ScalarSQTwoElectronOperator>> {

    // The type of transformation that is naturally associated to the second-quantized Hamiltonian.
    using Transformation = typename OperatorTraits<ScalarSQOneElectronOperator>::Transformation;

    // The type of the one-particle density matrix that is naturally associated to the second-quantized Hamiltonian.
    using OneDM = typename OperatorTraits<ScalarSQOneElectronOperator>::OneDM;

    // The type of the two-particle density matrix that is naturally associated to the second-quantized Hamiltonian.
    using TwoDM = typename OperatorTraits<ScalarSQOneElectronOperator>::TwoDM;
};


/*
 *  MARK: BasisTransformableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `BasisTransformable`.
 */
template <typename ScalarSQOneElectronOperator, typename ScalarSQTwoElectronOperator>
struct BasisTransformableTraits<SQHamiltonian<ScalarSQOneElectronOperator, ScalarSQTwoElectronOperator>> {

    // The type of the transformation for which the basis transformation should be defined.
    using Transformation = typename OperatorTraits<ScalarSQOneElectronOperator>::Transformation;
};


/*
 *  MARK: JacobiRotatableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `JacobiRotatable`.
 */
template <typename ScalarSQOneElectronOperator, typename ScalarSQTwoElectronOperator>
struct JacobiRotatableTraits<SQHamiltonian<ScalarSQOneElectronOperator, ScalarSQTwoElectronOperator>> {

    // The type of Jacobi rotation for which the Jacobi rotation should be defined.
    using JacobiRotationType = typename JacobiRotatableTraits<ScalarSQOneElectronOperator>::JacobiRotationType;
};


/*
 *  MARK: Convenience aliases
 */


// An `SQHamiltonian` related to restricted spin-orbitals. See `RestrictedSpinOrbitalTag`.
template <typename Scalar>
using RSQHamiltonian = SQHamiltonian<ScalarRSQOneElectronOperator<Scalar>, ScalarRSQTwoElectronOperator<Scalar>>;


// An `SQHamiltonian` related to unrestricted spin-orbitals. See `UnrestrictedSpinOrbitalTag`.
template <typename Scalar>
using USQHamiltonian = SQHamiltonian<ScalarUSQOneElectronOperator<Scalar>, ScalarUSQTwoElectronOperator<Scalar>>;


// An `SQHamiltonian` related to general spinors. See `GeneralSpinorTag`.
template <typename Scalar>
using GSQHamiltonian = SQHamiltonian<ScalarGSQOneElectronOperator<Scalar>, ScalarGSQTwoElectronOperator<Scalar>>;


}  // namespace GQCP
