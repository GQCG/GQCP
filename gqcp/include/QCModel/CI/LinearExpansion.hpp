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


#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Basis/SpinorBasis/GSpinorBasis.hpp"
#include "Basis/SpinorBasis/RSpinOrbitalBasis.hpp"
#include "Basis/SpinorBasis/USpinOrbitalBasis.hpp"
#include "Basis/Transformations/RTransformation.hpp"
#include "DensityMatrix/Orbital1DM.hpp"
#include "DensityMatrix/Orbital2DM.hpp"
#include "DensityMatrix/SpinResolved1DM.hpp"
#include "DensityMatrix/SpinResolved2DM.hpp"
#include "Mathematical/Representation/Matrix.hpp"
#include "ONVBasis/SpinResolvedONV.hpp"
#include "ONVBasis/SpinResolvedONVBasis.hpp"
#include "ONVBasis/SpinResolvedSelectedONVBasis.hpp"
#include "Utilities/aliases.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/range/adaptors.hpp>

#include <type_traits>


namespace GQCP {


/**
 *  A class that represents a linear expansion inside an ONV basis.
 * 
 *  @tparam _ONVBasis           The type of ONV basis.
 */
template <typename _ONVBasis>
class LinearExpansion {
public:
    // The type of the ONV basis.
    using ONVBasis = _ONVBasis;


private:
    // The ONV basis with respect to which the coefficients are defined
    ONVBasis onv_basis;

    // The expansion coefficients.
    VectorX<double> m_coefficients;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  Construct a linear expansion inside the given ONV basis, with corresponding expansion coefficients.
     *
     *  @param onv_basis            The ONV basis with respect to which the coefficients are defined
     *  @param coefficients         The expansion coefficients.
     */
    LinearExpansion(const ONVBasis& onv_basis, const VectorX<double>& coefficients) :
        onv_basis {onv_basis},
        m_coefficients {coefficients} {}


    /*
     *  The default constructor.
     */
    LinearExpansion() = default;


    /*
     *  MARK: Named constructors
     */

    /**
     *  Create a linear expansion with a normalized coefficient vector (i.e. all the coefficients are equal).
     * 
     *  @param onv_basis            The ONV basis with respect to which the coefficients are defined
     * 
     *  @return A constant LinearExpansion.
     */
    static LinearExpansion<ONVBasis> Constant(const ONVBasis& onv_basis) {

        VectorX<double> coefficients = VectorX<double>::Ones(onv_basis.dimension());
        coefficients.normalize();

        return LinearExpansion<ONVBasis>(onv_basis, coefficients);
    }


    /**
     *  Create a linear expansion that represents the Hartree-Fock wave function.
     * 
     *  @param onv_basis            The ONV basis with respect to which the coefficients are defined.
     * 
     *  @return a LinearExpansion
     */
    static LinearExpansion<ONVBasis> HartreeFock(const ONVBasis& onv_basis) {

        VectorX<double> coefficients = VectorX<double>::Unit(onv_basis.dimension(), 0);  // The first ONV in the ONV basis is expected to be the HF determinant.

        return LinearExpansion<ONVBasis>(onv_basis, coefficients);
    }


    /**
     *  Create a normalized linear expansion inside a given ONV basis with possibly non-normalized coefficients.
     * 
     *  @param onv_basis            The ONV basis with respect to which the coefficients are defined.
     *  @param coefficients         the expansion coefficients
     * 
     *  @return A normalized LinearExpansion.
     */
    static LinearExpansion<ONVBasis> Normalized(const ONVBasis& onv_basis, const VectorX<double>& coefficients) {

        // Normalize the coefficients if they aren't.
        return LinearExpansion<ONVBasis>(onv_basis,
                                         std::abs(coefficients.norm() - 1.0) > 1.0e-12 ? coefficients.normalized() : coefficients);
    }


    /**
     *  Create a linear expansion with a random, normalized coefficient vector, with coefficients uniformly distributed in [-1, +1] before any normalization.
     * 
     *  @param onv_basis            the ONV basis with respect to which the coefficients are defined
     * 
     *  @return A random LinearExpansion.
     */
    static LinearExpansion<ONVBasis> Random(const ONVBasis& onv_basis) {

        VectorX<double> coefficients = VectorX<double>::Random(onv_basis.dimension());
        coefficients.normalize();

        return LinearExpansion<ONVBasis>(onv_basis, coefficients);
    }


    /**
     *  Create a linear expansion by reading in a GAMESS-US file.
     * 
     *  @param GAMESSUS_filename      The name of the GAMESS-US file that contains the spin-resolved selected wave function expansion.
     * 
     *  @return The corresponding spin-resolved selected linear expansion from a given GAMESS-US file
     */
    template <typename Z = ONVBasis>
    static enable_if_t<std::is_same<Z, SpinResolvedSelectedONVBasis>::value, LinearExpansion<Z>> FromGAMESSUS(const std::string& GAMESSUS_filename) {

        // If the filename isn't properly converted into an input file stream, we assume the user supplied a wrong file.
        std::ifstream input_file_stream {GAMESSUS_filename};
        if (!input_file_stream.good()) {
            throw std::runtime_error("LinearExpansionReader(std::string): The provided GAMESS file is illegible. Maybe you specified a wrong path?");
        }

        std::ifstream input_file_stream_count {GAMESSUS_filename};  // made to count the expansion size


        // Do the actual parsing.

        // Read in dummy lines up until we actually get to the ONVs and coefficients.
        std::string line;
        std::string buffer;  // dummy for the counting stream
        while (std::getline(input_file_stream, line)) {
            std::getline(input_file_stream_count, buffer);


            if (line.find("ALPHA") != std::string::npos && line.find("BETA") != std::string::npos && line.find("COEFFICIENT") != std::string::npos) {  // if find() returns an index that's different from the 'not-found' index

                // This line should have dashes and we skip it
                std::getline(input_file_stream, line);
                std::getline(input_file_stream_count, buffer);
                break;
            }
        }

        size_t space_size = 0;
        // Count the number of configurations.
        while (std::getline(input_file_stream_count, buffer)) {
            if (buffer.length() != 0 | buffer[0] != '\n') {
                space_size++;
            }
        }

        VectorX<double> coefficients = VectorX<double>::Zero(space_size);


        std::getline(input_file_stream, line);

        // Read the first line containing the configurations.
        std::vector<std::string> splitted_line;
        boost::split(splitted_line, line, boost::is_any_of("|"));  // split on '|'


        // Create an alpha ONV for the first field
        std::string trimmed_alpha = boost::algorithm::trim_copy(splitted_line[0]);

        // Create a beta ONV for the second field
        std::string trimmed_beta = boost::algorithm::trim_copy(splitted_line[1]);


        size_t index_count = 0;  // counts the number of configurations added
        coefficients(index_count) = std::stod(splitted_line[2]);


        // Parse the trimmed ONV strings into boost::dynamic_bitset to use its functionality.
        std::string reversed_alpha {trimmed_alpha.rbegin(), trimmed_alpha.rend()};
        std::string reversed_beta {trimmed_beta.rbegin(), trimmed_beta.rend()};

        boost::dynamic_bitset<> alpha_transfer {reversed_alpha};
        boost::dynamic_bitset<> beta_transfer {reversed_beta};

        size_t K = alpha_transfer.size();
        size_t N_alpha = alpha_transfer.count();
        size_t N_beta = beta_transfer.count();

        SpinResolvedSelectedONVBasis onv_basis {K, N_alpha, N_beta};
        onv_basis.expandWith(SpinResolvedONV::FromString(reversed_alpha, reversed_beta));


        // Read in the ONVs and the coefficients by splitting the line on '|', and then trimming whitespace.
        // In the GAMESS-US format, the bit strings are ordered in reverse.
        while (std::getline(input_file_stream, line)) {

            index_count++;
            boost::split(splitted_line, line, boost::is_any_of("|"));  // split on '|'


            // Create an alpha SpinUnresolvedONV for the first field
            trimmed_alpha = boost::algorithm::trim_copy(splitted_line[0]);
            if (trimmed_alpha.length() != K) {
                throw std::invalid_argument("LinearExpansionReader(std::string): One of the provided alpha ONVs does not have the correct number of orbitals.");
            }
            reversed_alpha = std::string(trimmed_alpha.rbegin(), trimmed_alpha.rend());

            // Create a beta SpinUnresolvedONV for the second field
            trimmed_beta = boost::algorithm::trim_copy(splitted_line[1]);
            if (trimmed_beta.length() != K) {
                throw std::invalid_argument("LinearExpansionReader(std::string): One of the provided beta ONVs does not have the correct number of orbitals.");
            }
            reversed_beta = std::string(trimmed_beta.rbegin(), trimmed_beta.rend());


            // Create a double for the third field
            coefficients(index_count) = std::stod(splitted_line[2]);
            onv_basis.expandWith(SpinResolvedONV::FromString(reversed_alpha, reversed_beta));

        }  // while getline

        return LinearExpansion<SpinResolvedSelectedONVBasis>(onv_basis, coefficients);
    }


    /**
     *  Create the linear expansion of the given spin-resolved ONV that is expressed in the given USpinOrbitalBasis, by projection onto the spin-resolved ONVs expressed with respect to the given RSpinOrbitalBasis.
     * 
     *  @param onv                      A spin-resolved ONV expressed with respect to an unrestricted spin-orbital basis.
     *  @param r_spinor_basis           The restricted spin-orbital basis that is used to define the resulting linear expansion of ONVs against.
     *  @param u_spinor_basis           The unrestricted spin-orbital basis against which the given ONV is expressed.
     * 
     *  @return A linear expansion inside a spin-resolved ONV basis.
     */
    template <typename Z = ONVBasis>
    static enable_if_t<std::is_same<Z, SpinResolvedONVBasis>::value, LinearExpansion<Z>> FromONVProjection(const SpinResolvedONV& onv, const RSpinOrbitalBasis<double, GTOShell>& r_spinor_basis, const USpinOrbitalBasis<double, GTOShell>& u_spinor_basis) {

        // Determine the overlap matrices of the underlying scalar orbital bases, which is needed later on.
        auto S_r = r_spinor_basis.overlap();                  // the overlap matrix of the restricted MOs/spin-orbitals
        S_r.transform(r_spinor_basis.expansion().inverse());  // now in AO basis

        auto S_u = u_spinor_basis.overlap();                  // The overlap matrix of the unrestricted spin-orbitals.
        S_u.transform(u_spinor_basis.expansion().inverse());  // Now in AO basis.

        if (!(S_r.parameters().isApprox(S_u.alpha().parameters(), 1.0e-08)) || !(S_r.parameters().isApprox(S_u.beta().parameters(), 1.0e-08))) {
            throw std::invalid_argument("LinearExpansion::FromONVProjection(const SpinResolvedONV&, const RSpinOrbitalBasis<double, GTOShell>&, const USpinOrbitalBasis<double, GTOShell>&): The given spinor bases are not expressed using the same scalar orbital basis.");
        }


        // Prepare some parameters.
        const auto& C_restricted = r_spinor_basis.expansion();

        const auto C_unrestricted = u_spinor_basis.expansion();


        // Set up the required spin-resolved ONV basis.
        const auto K = onv.numberOfSpatialOrbitals(Spin::alpha);  // assume equal numbers of spin-orbitals for alpha- and beta-electrons
        const auto N_alpha = onv.numberOfElectrons(Spin::alpha);
        const auto N_beta = onv.numberOfElectrons(Spin::beta);
        const SpinResolvedONVBasis onv_basis {K, N_alpha, N_beta};


        // Determine the coefficients through calculating the overlap between two ONVs.
        VectorX<double> coefficients = VectorX<double>::Zero(onv_basis.dimension());

        onv_basis.forEach([&onv, &C_unrestricted, &C_restricted, &S_r, &coefficients, &onv_basis](const SpinUnresolvedONV& alpha_onv, const size_t I_alpha, const SpinUnresolvedONV& beta_onv, const size_t I_beta) {
            const SpinResolvedONV onv_on {alpha_onv, beta_onv};  // the spin-resolved ONV that should be projected 'on'

            const auto coefficient = onv.calculateProjection(onv_on, C_unrestricted, C_restricted, S_r.parameters());
            const auto address = onv_basis.compoundAddress(I_alpha, I_beta);

            coefficients(address) = coefficient;
        });

        return LinearExpansion<Z>(onv_basis, coefficients);
    }


    /**
     *  Create the linear expansion of the given spin-unresolved ONV that is expressed in the given GSpinorBasis, by projection onto the spin-resolved ONVs expressed with respect to another given GSpinorBasis.
     * 
     *  @param onv_of                   A spin-unresolved ONV expressed with respect to a general spinor basis.
     *  @param spinor_basis_on          The general spinor basis that is used to define the resulting linear expansion of ONVs against.
     *  @param spinor_basis_of          The general spinor basis against which the given ONV is expressed.
     * 
     *  @return A linear expansion inside a spin-unresolved ONV basis.
     */
    template <typename Z = ONVBasis>
    static enable_if_t<std::is_same<Z, SpinUnresolvedONVBasis>::value, LinearExpansion<Z>> FromONVProjection(const SpinUnresolvedONV& onv_of, const GSpinorBasis<double, GTOShell>& spinor_basis_on, const GSpinorBasis<double, GTOShell>& spinor_basis_of) {

        // Determine the overlap matrices of the underlying scalar orbital bases, which is needed later on.
        auto S_on = spinor_basis_on.overlap();
        S_on.transform(spinor_basis_on.expansion().inverse());  // now in AO basis

        auto S_of = spinor_basis_of.overlap();
        S_of.transform(spinor_basis_of.expansion().inverse());  // now in AO basis

        if (!(S_on.parameters().isApprox(S_of.parameters(), 1.0e-08))) {
            throw std::invalid_argument("LinearExpansion::FromONVProjection(const SpinUnresolvedONV&, const RSpinOrbitalBasis<double, GTOShell>&, const GSpinorBasis<double, GTOShell>&): The given spinor bases are not expressed using the same scalar orbital basis.");
        }


        // Prepare some parameters.
        const auto& C_on = spinor_basis_on.expansion();
        const auto& C_of = spinor_basis_of.expansion();


        // Set up the required spin-resolved ONV basis.
        const auto M = onv_of.numberOfSpinors();
        const auto N = onv_of.numberOfElectrons();
        const SpinUnresolvedONVBasis onv_basis {M, N};


        // Determine the coefficients through calculating the overlap between two ONVs.
        VectorX<double> coefficients = VectorX<double>::Zero(onv_basis.dimension());
        onv_basis.forEach([&onv_of, &C_on, &C_of, &S_on, &coefficients](const SpinUnresolvedONV& onv_on, const size_t I) {
            coefficients(I) = onv_of.calculateProjection(onv_on, C_of, C_on, S_on.parameters());
        });

        return LinearExpansion<Z>(onv_basis, coefficients);
    }


    /*
     *  MARK: Access
     */

    /**
     *  Access a coefficient of the linear expansion.
     * 
     *  @param i    The index (address) of the coefficient that should be obtained.
     * 
     *  @return The i-th expansion coefficient.
     */
    double coefficient(const size_t i) const { return this->m_coefficients(i); }

    /**
     *  @return The expansion coefficients of this linear expansion wave function model.
     */
    const VectorX<double>& coefficients() const { return this->m_coefficients; }

    /**
     *  @return The ONV basis that is related to this linear expansion wave function model.
     */
    const ONVBasis& onvBasis() const { return onv_basis; }


    /*
     *  MARK: Basis transformations
     */

    /**
     *  Update the expansion coefficients of this linear expansion so that they correspond to the situation after a transformation of the underlying spinor basis with the given basis transformation.
     *
     *  @param T            The transformation between the old and the new restricted spin-orbital basis.
     * 
     *  @note This method is only available for the full spin-resolved ONV basis.
     *  @note This algorithm was implemented from a description in Helgaker2000.
     */
    template <typename Z = ONVBasis>
    enable_if_t<std::is_same<Z, SpinResolvedONVBasis>::value> basisTransform(const RTransformation<double>& T) {

        const auto K = onv_basis.numberOfOrbitals();  // number of spatial orbitals
        if (K != T.numberOfOrbitals()) {
            throw std::invalid_argument("LinearExpansion::basisTransform(const RTransformation<double>&): The number of spatial orbitals does not match the dimension of the transformation matrix.");
        }


        // LU-decompose the transformation matrix LU decomposition for T
        const auto& lu_decomposition = T.matrix().noPivotLUDecompose();

        SquareMatrix<double> L = SquareMatrix<double>::Identity(K);
        L.triangularView<Eigen::StrictlyLower>() = lu_decomposition[0];

        SquareMatrix<double> U = SquareMatrix<double>(lu_decomposition[1].triangularView<Eigen::Upper>());
        SquareMatrix<double> U_inv = U.inverse();


        // Calculate t (the operator which allows per-orbital transformation of the wave function)
        SquareMatrix<double> t = SquareMatrix<double>::Identity(K) - L + U_inv;


        // Set up spin-unresolved ONV basis variables for the loops over the ONVs
        const SpinUnresolvedONVBasis& alpha_onv_basis = onv_basis.alpha();
        const SpinUnresolvedONVBasis& beta_onv_basis = onv_basis.beta();

        auto dim_alpha = alpha_onv_basis.dimension();
        auto dim_beta = beta_onv_basis.dimension();
        auto N_alpha = alpha_onv_basis.numberOfElectrons();
        auto N_beta = beta_onv_basis.numberOfElectrons();


        /** 
         *  The transformation of the expansion coefficients is adapted from Helgaker2000, chapter 11.9.
         *  For every orbital, a set of correction coefficients will be calculated (Delta C in Helgaker), to update the current coefficients.
         */

        VectorX<double> current_coefficients = this->m_coefficients;  // coefficients will be updated after each orbital transform (C^(n-1)) in Helgaker
        VectorX<double> correction_coefficients = VectorX<double>::Zero(onv_basis.dimension());


        for (size_t m = 0; m < K; m++) {  // iterate over all orbitals

            // Perform alpha and beta CI iterations.

            // 1) Alpha-branch
            SpinUnresolvedONV alpha = alpha_onv_basis.constructONVFromAddress(0);
            for (size_t I_alpha = 0; I_alpha < dim_alpha; I_alpha++) {
                if (!alpha.isOccupied(m)) {
                    for (size_t e1 = 0; e1 < N_alpha; e1++) {    // e1 (electron 1) loops over the (number of) electrons
                        size_t p = alpha.occupationIndexOf(e1);  // retrieve the index of a given electron

                        if (p < m) {
                            size_t address = I_alpha - alpha_onv_basis.vertexWeight(p, e1 + 1);
                            size_t e2 = e1 + 1;
                            size_t q = p + 1;
                            int sign = 1;

                            alpha_onv_basis.shiftUntilNextUnoccupiedOrbital<1>(alpha, address, q, e2, sign);
                            while (q != m) {
                                q++;
                                alpha_onv_basis.shiftUntilNextUnoccupiedOrbital<1>(alpha, address, q, e2, sign);
                            }

                            address += alpha_onv_basis.vertexWeight(q, e2);

                            for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {
                                correction_coefficients(I_alpha * dim_beta + I_beta) += sign * t(p, m) * current_coefficients(address * dim_beta + I_beta);
                            }
                        }

                        if (p > m) {
                            size_t address = I_alpha - alpha_onv_basis.vertexWeight(p, e1 + 1);
                            size_t e2 = e1 - 1;
                            size_t q = p - 1;
                            int sign = 1;

                            alpha_onv_basis.shiftUntilPreviousUnoccupiedOrbital<1>(alpha, address, q, e2, sign);
                            while (q != m) {
                                q--;
                                alpha_onv_basis.shiftUntilPreviousUnoccupiedOrbital<1>(alpha, address, q, e2, sign);
                            }

                            address += alpha_onv_basis.vertexWeight(q, e2 + 2);
                            for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {
                                correction_coefficients(I_alpha * dim_beta + I_beta) += sign * t(p, m) * current_coefficients(address * dim_beta + I_beta);
                            }
                        }
                    }

                } else {  // if orbital m is occupied we can perform an in-place operation
                    for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {
                        correction_coefficients(I_alpha * dim_beta + I_beta) += (t(m, m) - 1) * current_coefficients(I_alpha * dim_beta + I_beta);
                    }
                }


                if (I_alpha < dim_alpha - 1) {  // prevent the last permutation from occurring
                    alpha_onv_basis.transformONVToNextPermutation(alpha);
                }
            }

            current_coefficients += correction_coefficients;
            correction_coefficients.setZero();

            // 2) Beta-branch
            SpinUnresolvedONV beta = beta_onv_basis.constructONVFromAddress(0);

            for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {
                if (!beta.isOccupied(m)) {
                    for (size_t e1 = 0; e1 < N_beta; e1++) {    // e1 (electron 1) loops over the (number of) electrons
                        size_t p = beta.occupationIndexOf(e1);  // retrieve the index of a given electron

                        if (p < m) {
                            size_t address = I_beta - beta_onv_basis.vertexWeight(p, e1 + 1);
                            size_t e2 = e1 + 1;
                            size_t q = p + 1;
                            int sign = 1;

                            beta_onv_basis.shiftUntilNextUnoccupiedOrbital<1>(beta, address, q, e2, sign);
                            while (q != m) {
                                q++;
                                beta_onv_basis.shiftUntilNextUnoccupiedOrbital<1>(beta, address, q, e2, sign);
                            }

                            address += beta_onv_basis.vertexWeight(q, e2);

                            for (size_t I_alpha = 0; I_alpha < dim_alpha; I_alpha++) {
                                correction_coefficients(I_alpha * dim_beta + I_beta) += sign * t(p, m) * current_coefficients(I_alpha * dim_beta + address);
                            }
                        }

                        if (p > m) {
                            size_t address = I_beta - beta_onv_basis.vertexWeight(p, e1 + 1);
                            size_t e2 = e1 - 1;
                            size_t q = p - 1;

                            int sign = 1;

                            beta_onv_basis.shiftUntilPreviousUnoccupiedOrbital<1>(beta, address, q, e2, sign);
                            while (q != m) {
                                q--;
                                beta_onv_basis.shiftUntilPreviousUnoccupiedOrbital<1>(beta, address, q, e2, sign);
                            }

                            address += beta_onv_basis.vertexWeight(q, e2 + 2);

                            for (size_t I_alpha = 0; I_alpha < dim_alpha; I_alpha++) {
                                correction_coefficients(I_alpha * dim_beta + I_beta) += sign * t(p, m) * current_coefficients(I_alpha * dim_beta + address);
                            }
                        }
                    }

                } else {  // if orbital m is occupied we can perform an in-place operation
                    for (size_t I_alpha = 0; I_alpha < dim_alpha; I_alpha++) {
                        correction_coefficients(I_alpha * dim_beta + I_beta) += (t(m, m) - 1) * current_coefficients(I_alpha * dim_beta + I_beta);
                    }
                }

                if (I_beta < dim_beta - 1) {  // prevent the last permutation from occurring
                    beta_onv_basis.transformONVToNextPermutation(beta);
                }
            }

            current_coefficients += correction_coefficients;
            correction_coefficients.setZero();
        }
        this->m_coefficients = current_coefficients;
    }


    /*
     *  MARK: Density matrices for spin-unresolved ONV bases
     */

    /*
     *  Calculate general one-electron density matrix for a spin-unresolved wave function expansion.
     * 
     *  @return The generalized one-electron density matrix.
     */
    template <typename Z = ONVBasis>
    enable_if_t<std::is_same<Z, SpinUnresolvedONVBasis>::value, G1DM<double>> calculate1DM() const {

        // Prepare some variables.
        const auto M = this->onv_basis.numberOfOrbitals();
        const auto N = this->onv_basis.numberOfElectrons();
        const auto dim = onv_basis.dimension();  // Dimension of the SpinUnresolvedONVBasis = number of SpinUnresolvedONVs.

        GQCP::G1DM<double> D = GQCP::G1DM<double>::Zero(M);

        SpinUnresolvedONV onv = onv_basis.constructONVFromAddress(0);  // Start with ONV with address 0.
        for (size_t J = 0; J < dim; J++) {                             // Loops over all possible ONV indices.

            // Loop over electrons that can be annihilated in an ONV.
            for (size_t e1 = 0; e1 < N; e1++) {

                // Create an ONVPath for each new ONV.
                ONVPath<SpinUnresolvedONVBasis> onv_path {onv_basis, onv};
                const auto c_J = this->coefficient(J);

                // Figure out the orbital index of the electron that will be annihilated.
                const auto q = onv.occupationIndexOf(e1);


                // The diagonal values are a result of annihilation-creation on the same orbital index and are thus the same as the initial ONV.
                D(q, q) += c_J * c_J;

                // For the non-diagonal values, we will create all possible matrix elements of the density matrix in the routine below.
                onv_path.annihilate(q, e1);

                // Stop the loop if 1) the path is finished, meaning that orbital index p is at M (the total number of orbitals) and 2) if the orbital index is out of bounds after left translation of a vertical arc.
                while (!onv_path.isFinished() && onv_path.isOrbitalIndexValid()) {

                    // Find the next unoccupied orbital, i.e. the next vertical arc in the path.
                    onv_path.leftTranslateDiagonalArcUntilVerticalArc();

                    // Calculate the address of the path if we would close it right now.
                    const auto I = onv_path.addressAfterCreation();
                    const auto c_I = this->coefficient(I);


                    // Add the density matrix elements.
                    const auto p = onv_path.orbitalIndex();

                    const double value = c_I * c_J;

                    D(p, q) += onv_path.sign() * value;
                    D(q, p) += onv_path.sign() * value;

                    // Move orbital index such that other unoccupied orbitals can be found within the loop.
                    onv_path.leftTranslateVerticalArc();
                }
            }

            // Prevent last ONV since there is no possibility for an electron to be annihilated anymore.
            if (J < dim - 1) {
                onv_basis.transformONVToNextPermutation(onv);
            }
        }

        return D;
    }


    /**
     *  Calculate an element of the N-electron density matrix.
     * 
     *  @param bra_indices      The indices of the orbitals that should be annihilated on the left (on the bra).
     *  @param ket_indices      The indices of the orbitals that should be annihilated on the right (on the ket).
     *
     *  @return An element of the N-DM, as specified by the given bra and ket indices. `calculateNDMElement({0, 1}, {2, 1})` would calculate an element of the 2-NDM d^{(2)} (0, 1, 1, 2) corresponding the operator string: `a^\dagger_0 a^\dagger_1 a_2 a_1`.
     * 
     *  @note This method is only enabled for linear expansions related to spin-unresolved ONV bases.
     */
    template <typename Z = ONVBasis>
    enable_if_t<std::is_same<Z, SpinUnresolvedONVBasis>::value, double> calculateNDMElement(const std::vector<size_t>& bra_indices, const std::vector<size_t>& ket_indices) const {

        // The ket indices should be reversed because the annihilators on the ket should be applied from right to left
        std::vector<size_t> ket_indices_reversed = ket_indices;
        std::reverse(ket_indices_reversed.begin(), ket_indices_reversed.end());


        double value = 0.0;
        int sign = 1;
        const size_t dim = this->onv_basis.dimension();


        SpinUnresolvedONV bra = this->onv_basis.constructONVFromAddress(0);
        size_t I = 0;
        while (I < dim) {  // loop over all bra addresses

            // Annihilate the bra on the bra indices
            if (!bra.annihilateAll(bra_indices, sign)) {  // if we can't annihilate, the bra doesn't change

                // Go to the beginning of the outer while loop with the next bra
                if (I < dim - 1) {  // prevent the last permutation from occurring
                    this->onv_basis.transformONVToNextPermutation(bra);
                    I++;
                    sign = 1;
                    continue;
                } else {
                    break;  // we have to jump out if we have looped over the whole bra dimension
                }
            }


            SpinUnresolvedONV ket = this->onv_basis.constructONVFromAddress(0);
            size_t J = 0;
            while (J < dim) {  // loop over all ket indices

                // Annihilate the ket on the ket indices
                if (!ket.annihilateAll(ket_indices_reversed, sign)) {  // if we can't annihilate, the ket doesn't change
                    // Go to the beginning of this (the inner) while loop with the next bra
                    if (J < dim - 1) {  // prevent the last permutation from occurring
                        this->onv_basis.transformONVToNextPermutation(ket);
                        J++;
                        sign = 1;
                        continue;
                    } else {
                        break;  // we have to jump out if we have looped over the whole ket dimension
                    }
                }

                if (bra == ket) {
                    value += sign * this->coefficient(I) * this->coefficient(J);
                }

                // Reset the previous ket annihilations and move to the next ket
                if (J == dim - 1) {  // prevent the last permutation from occurring
                    break;           // out of the J-loop
                }
                ket.createAll(ket_indices_reversed);
                this->onv_basis.transformONVToNextPermutation(ket);
                sign = 1;
                J++;
            }  // while J loop

            // Reset the previous bra annihilations and move to the next bra
            if (I == dim - 1) {  // prevent the last permutation from occurring
                break;           // out of the I-loop
            }
            bra.createAll(bra_indices);
            this->onv_basis.transformONVToNextPermutation(bra);
            sign = 1;
            I++;
        }  // while I loop

        return value;
    }


    /*
     *  MARK: Density matrices for spin-resolved ONV bases
     */

    /**
     *  Calculate the spin-resolved one-electron density matrix for a full spin-resolved wave function expansion.
     * 
     *  @return The spin-resolved 1-DM.
     */
    template <typename Z = ONVBasis>
    enable_if_t<std::is_same<Z, SpinResolvedONVBasis>::value, SpinResolved1DM<double>> calculateSpinResolved1DM() const {

        // Initialize as zero matrices
        size_t K = this->onv_basis.alpha().numberOfOrbitals();

        SpinResolved1DMComponent<double> D_aa = SpinResolved1DMComponent<double>::Zero(K);
        SpinResolved1DMComponent<double> D_bb = SpinResolved1DMComponent<double>::Zero(K);

        SpinUnresolvedONVBasis onv_basis_alpha = onv_basis.alpha();
        SpinUnresolvedONVBasis onv_basis_beta = onv_basis.beta();

        auto dim_alpha = onv_basis_alpha.dimension();
        auto dim_beta = onv_basis_beta.dimension();

        // ALPHA
        SpinUnresolvedONV spin_string_alpha = onv_basis_alpha.constructONVFromAddress(0);  // alpha spin string with address 0
        for (size_t I_alpha = 0; I_alpha < dim_alpha; I_alpha++) {                         // I_alpha loops over all the addresses of the alpha spin strings
            for (size_t p = 0; p < K; p++) {                                               // p loops over SOs
                int sign_p = 1;
                if (spin_string_alpha.annihilate(p, sign_p)) {  // if p is in I_alpha
                    double diagonal_contribution = 0;

                    // Diagonal contributions for the 1-DM, i.e. D_pp
                    // We are storing the alpha addresses as 'major', i.e. the total address I_alpha I_beta = I_alpha * dim_beta + I_beta
                    for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {
                        double c_I_alpha_I_beta = this->coefficient(I_alpha * dim_beta + I_beta);
                        diagonal_contribution += std::pow(c_I_alpha_I_beta, 2);
                    }
                    D_aa(p, p) += diagonal_contribution;

                    // Off-diagonal contributions for the 1-DM, i.e. D_pq (p!=q)
                    for (size_t q = 0; q < p; q++) {  // q < p loops over SOs
                        int sign_pq = sign_p;
                        if (spin_string_alpha.create(q, sign_pq)) {                         // if q is not occupied in I_alpha
                            size_t J_alpha = onv_basis_alpha.addressOf(spin_string_alpha);  // find all strings J_alpha that couple to I_alpha

                            double off_diagonal_contribution = 0;
                            for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {
                                double c_I_alpha_I_beta = this->coefficient(I_alpha * dim_beta + I_beta);  // alpha addresses are 'major'
                                double c_J_alpha_I_beta = this->coefficient(J_alpha * dim_beta + I_beta);
                                off_diagonal_contribution += c_I_alpha_I_beta * c_J_alpha_I_beta;
                            }

                            D_aa(p, q) += sign_pq * off_diagonal_contribution;
                            D_aa(q, p) += sign_pq * off_diagonal_contribution;  // add the symmetric contribution because we are looping over q < p

                            spin_string_alpha.annihilate(q);  // undo the previous creation

                        }  // create on q
                    }      // q loop

                    spin_string_alpha.create(p);  // undo the previous annihilation

                }  // annihilate on p
            }      // p loop

            if (I_alpha < dim_alpha - 1) {  // prevent the last permutation from occurring
                onv_basis_alpha.transformONVToNextPermutation(spin_string_alpha);
            }

        }  // I_alpha loop


        // BETA
        SpinUnresolvedONV spin_string_beta = onv_basis_beta.constructONVFromAddress(0);  // spin string with address 0
        for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {                           // I_beta loops over all the addresses of the spin strings
            for (size_t p = 0; p < K; p++) {                                             // p loops over SOs
                int sign_p = 1;
                if (spin_string_beta.annihilate(p, sign_p)) {  // if p is in I_beta
                    double diagonal_contribution = 0;

                    // Diagonal contributions for the 1-DM, i.e. D_pp
                    // We are storing the alpha addresses as 'major', i.e. the total address I_alpha I_beta = I_alpha * dim_beta + I_beta
                    for (size_t I_alpha = 0; I_alpha < dim_alpha; I_alpha++) {
                        double c_I_alpha_I_beta = this->coefficient(I_alpha * dim_beta + I_beta);
                        diagonal_contribution += std::pow(c_I_alpha_I_beta, 2);
                    }

                    D_bb(p, p) += diagonal_contribution;

                    // Off-diagonal contributions for the 1-DM
                    for (size_t q = 0; q < p; q++) {  // q < p loops over SOs
                        int sign_pq = sign_p;
                        if (spin_string_beta.create(q, sign_pq)) {                       // if q is not in I_beta
                            size_t J_beta = onv_basis_beta.addressOf(spin_string_beta);  // find all strings J_beta that couple to I_beta

                            double off_diagonal_contribution = 0;
                            for (size_t I_alpha = 0; I_alpha < dim_alpha; I_alpha++) {
                                double c_I_alpha_I_beta = this->coefficient(I_alpha * dim_beta + I_beta);  // alpha addresses are 'major'
                                double c_I_alpha_J_beta = this->coefficient(I_alpha * dim_beta + J_beta);
                                off_diagonal_contribution += c_I_alpha_I_beta * c_I_alpha_J_beta;
                            }
                            D_bb(p, q) += sign_pq * off_diagonal_contribution;
                            D_bb(q, p) += sign_pq * off_diagonal_contribution;  // add the symmetric contribution because we are looping over q < p

                            spin_string_beta.annihilate(q);  // undo the previous creation

                        }  // create on q
                    }      // loop over q

                    spin_string_beta.create(p);  // undo the previous annihilation

                }  // annihilate on p
            }      // loop over p

            if (I_beta < dim_beta - 1) {  // prevent the last permutation from occurring
                onv_basis_beta.transformONVToNextPermutation(spin_string_beta);
            }

        }  // I_beta loop
        return SpinResolved1DM<double>(D_aa, D_bb);
    }


    /**
     *  Calculate the spin-resolved two-electron density matrix for a full spin-resolved wave function expansion.
     * 
     *  @return The spin-resolved 2-DM.
     */
    template <typename Z = ONVBasis>
    enable_if_t<std::is_same<Z, SpinResolvedONVBasis>::value, SpinResolved2DM<double>> calculateSpinResolved2DM() const {

        // KISS implementation of the 2-DMs (no symmetry relations are used yet)

        SpinUnresolvedONVBasis onv_basis_alpha = onv_basis.alpha();
        SpinUnresolvedONVBasis onv_basis_beta = onv_basis.beta();

        auto dim_alpha = onv_basis_alpha.dimension();
        auto dim_beta = onv_basis_beta.dimension();

        // Initialize as zero matrices
        size_t K = this->onv_basis.alpha().numberOfOrbitals();

        PureSpinResolved2DMComponent<double> d_aaaa = PureSpinResolved2DMComponent<double>::Zero(K);
        MixedSpinResolved2DMComponent<double> d_aabb = MixedSpinResolved2DMComponent<double>::Zero(K);
        PureSpinResolved2DMComponent<double> d_bbbb = PureSpinResolved2DMComponent<double>::Zero(K);

        // ALPHA-ALPHA-ALPHA-ALPHA
        SpinUnresolvedONV spin_string_alpha_aaaa = onv_basis_alpha.constructONVFromAddress(0);  // spin string with address 0
        for (size_t I_alpha = 0; I_alpha < dim_alpha; I_alpha++) {                              // I_alpha loops over all the addresses of the alpha spin strings

            for (size_t p = 0; p < K; p++) {  // p loops over SOs
                int sign_p = 1;               // sign of the operator a_p

                if (spin_string_alpha_aaaa.annihilate(p, sign_p)) {  // if p is not in I_alpha

                    for (size_t r = 0; r < K; r++) {  // r loops over SOs
                        int sign_pr = sign_p;         // sign of the operator a_r a_p

                        if (spin_string_alpha_aaaa.annihilate(r, sign_pr)) {  // if r is not in I_alpha

                            for (size_t s = 0; s < K; s++) {  // s loops over SOs
                                int sign_prs = sign_pr;       // sign of the operator a^dagger_s a_r a_p

                                if (spin_string_alpha_aaaa.create(s, sign_prs)) {  // if s is in I_alpha

                                    for (size_t q = 0; q < K; q++) {  // q loops over SOs
                                        int sign_prsq = sign_prs;     // sign of the operator a^dagger_q a^dagger_s a_r a_p

                                        if (spin_string_alpha_aaaa.create(q, sign_prsq)) {                       // if q is not in I_alpha
                                            size_t J_alpha = onv_basis_alpha.addressOf(spin_string_alpha_aaaa);  // address of the coupling string

                                            double contribution = 0.0;
                                            for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {
                                                double c_I_alpha_I_beta = this->coefficient(I_alpha * dim_beta + I_beta);  // alpha addresses are 'major'
                                                double c_J_alpha_I_beta = this->coefficient(J_alpha * dim_beta + I_beta);
                                                contribution += c_I_alpha_I_beta * c_J_alpha_I_beta;
                                            }


                                            d_aaaa(p, q, r, s) += sign_prsq * contribution;

                                            spin_string_alpha_aaaa.annihilate(q);  // undo the previous creation
                                        }
                                    }  // loop over q

                                    spin_string_alpha_aaaa.annihilate(s);  // undo the previous creation
                                }
                            }  // loop over s

                            spin_string_alpha_aaaa.create(r);  // undo the previous annihilation
                        }
                    }  // loop over r

                    spin_string_alpha_aaaa.create(p);  // undo the previous annihilation
                }
            }  // loop over p

            if (I_alpha < dim_alpha - 1) {  // prevent the last permutation from occurring
                onv_basis_alpha.transformONVToNextPermutation(spin_string_alpha_aaaa);
            }

        }  // loop over I_alpha


        // ALPHA-ALPHA-BETA-BETA
        SpinUnresolvedONV spin_string_alpha_aabb = onv_basis_alpha.constructONVFromAddress(0);  // spin string with address 0
        for (size_t I_alpha = 0; I_alpha < dim_alpha; I_alpha++) {                              // I_alpha loops over all the addresses of the alpha spin strings

            for (size_t p = 0; p < K; p++) {  // p loops over SOs
                int sign_p = 1;               // sign of the operator a_p_alpha

                if (spin_string_alpha_aabb.annihilate(p, sign_p)) {  // if p is in I_alpha

                    for (size_t q = 0; q < K; q++) {  // q loops over SOs
                        int sign_pq = sign_p;         // sign of the operator a^dagger_p_alpha a_p_alpha

                        if (spin_string_alpha_aabb.create(q, sign_pq)) {                         // if q is not in I_alpha
                            size_t J_alpha = onv_basis_alpha.addressOf(spin_string_alpha_aabb);  // the string that couples to I_alpha


                            SpinUnresolvedONV spin_string_beta_aabb = onv_basis_beta.constructONVFromAddress(0);  // spin string with address 0
                            for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {                                // I_beta loops over all addresses of beta spin strings

                                for (size_t r = 0; r < K; r++) {  // r loops over all SOs
                                    int sign_r = 1;               // sign of the operator a_r_beta

                                    if (spin_string_beta_aabb.annihilate(r, sign_r)) {

                                        for (size_t s = 0; s < K; s++) {  // s loops over all SOs
                                            int sign_rs = sign_r;         // sign of the operator a^dagger_s_beta a_r_beta

                                            if (spin_string_beta_aabb.create(s, sign_rs)) {
                                                size_t J_beta = onv_basis_beta.addressOf(spin_string_beta_aabb);  // the string that couples to I_beta

                                                double c_I_alpha_I_beta = this->coefficient(I_alpha * dim_beta + I_beta);  // alpha addresses are 'major'
                                                double c_J_alpha_J_beta = this->coefficient(J_alpha * dim_beta + J_beta);
                                                d_aabb(p, q, r, s) += sign_pq * sign_rs * c_I_alpha_I_beta * c_J_alpha_J_beta;

                                                spin_string_beta_aabb.annihilate(s);  // undo the previous creation
                                            }
                                        }  // loop over s


                                        spin_string_beta_aabb.create(r);  // undo the previous annihilation
                                    }

                                }  // loop over r

                                if (I_beta < dim_beta - 1) {  // prevent the last permutation from occurring
                                    onv_basis_beta.transformONVToNextPermutation(spin_string_beta_aabb);
                                }

                            }  // loop over beta addresses

                            spin_string_alpha_aabb.annihilate(q);  // undo the previous creation
                        }
                    }  // loop over q

                    spin_string_alpha_aabb.create(p);  // undo the previous annihilation
                }
            }  // loop over p

            if (I_alpha < dim_alpha - 1) {  // prevent the last permutation from occurring
                onv_basis_alpha.transformONVToNextPermutation(spin_string_alpha_aabb);
            }

        }  // loop over alpha addresses


        // BETA-BETA-ALPHA-ALPHA
        // We know that d^aabb_pqrs = d^bbaa_rspq
        Eigen::array<int, 4> shuffle {2, 3, 0, 1};  // array specifying the axes that should be swapped
        MixedSpinResolved2DMComponent<double> d_bbaa {d_aabb.Eigen().shuffle(shuffle)};


        // BETA-BETA-BETA-BETA
        SpinUnresolvedONV spin_string_beta_bbbb = onv_basis_beta.constructONVFromAddress(0);  // spin string with address 0
        for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {                                // I_beta loops over all the addresses of the beta spin strings

            for (size_t p = 0; p < K; p++) {  // p loops over SOs
                int sign_p = 1;               // sign of the operator a_p

                if (spin_string_beta_bbbb.annihilate(p, sign_p)) {  // if p is not in I_beta

                    for (size_t r = 0; r < K; r++) {  // r loops over SOs
                        int sign_pr = sign_p;         // sign of the operator a_r a_p

                        if (spin_string_beta_bbbb.annihilate(r, sign_pr)) {  // if r is not in I_beta

                            for (size_t s = 0; s < K; s++) {  // s loops over SOs
                                int sign_prs = sign_pr;       // sign of the operator a^dagger_s a_r a_p

                                if (spin_string_beta_bbbb.create(s, sign_prs)) {  // if s is in I_beta

                                    for (size_t q = 0; q < K; q++) {  // q loops over SOs
                                        int sign_prsq = sign_prs;     // sign of the operator a^dagger_q a^dagger_s a_r a_p

                                        if (spin_string_beta_bbbb.create(q, sign_prsq)) {                     // if q is not in I_beta
                                            size_t J_beta = onv_basis_beta.addressOf(spin_string_beta_bbbb);  // address of the coupling string

                                            double contribution = 0.0;
                                            for (size_t I_alpha = 0; I_alpha < dim_alpha; I_alpha++) {
                                                double c_I_alpha_I_beta = this->coefficient(I_alpha * dim_beta + I_beta);  // alpha addresses are 'major'
                                                double c_I_alpha_J_beta = this->coefficient(I_alpha * dim_beta + J_beta);
                                                contribution += c_I_alpha_I_beta * c_I_alpha_J_beta;
                                            }


                                            d_bbbb(p, q, r, s) += sign_prsq * contribution;

                                            spin_string_beta_bbbb.annihilate(q);  // undo the previous creation
                                        }
                                    }  // loop over q

                                    spin_string_beta_bbbb.annihilate(s);  // undo the previous creation
                                }
                            }  // loop over s

                            spin_string_beta_bbbb.create(r);  // undo the previous annihilation
                        }
                    }  // loop over r

                    spin_string_beta_bbbb.create(p);  // undo the previous annihilation
                }
            }  // loop over p

            if (I_beta < dim_beta - 1) {  // prevent the last permutation from occurring
                onv_basis_beta.transformONVToNextPermutation(spin_string_beta_bbbb);
            }

        }  // loop over I_beta

        return SpinResolved2DM<double> {d_aaaa, d_aabb, d_bbaa, d_bbbb};
    }


    /**
     *  Calculate the one-electron density matrix for a full spin-resolved wave function expansion.
     * 
     *  @return The orbital (total, spin-summed) 1-DM
     */
    template <typename Z = ONVBasis>
    enable_if_t<std::is_same<Z, SpinResolvedONVBasis>::value, Orbital1DM<double>> calculate1DM() const { return this->calculateSpinResolved1DM().orbitalDensity(); }


    /**
     *  Calculate the two-electron density matrix for a full spin-resolved wave function expansion.
     * 
     *  @return The orbital (total, spin-summed) 2-DM.
     */
    template <typename Z = ONVBasis>
    enable_if_t<std::is_same<Z, SpinResolvedONVBasis>::value, Orbital2DM<double>> calculate2DM() const { return this->calculateSpinResolved2DM().orbitalDensity(); }


    /*
     *  MARK: Density matrices for seniority-zero ONV bases
     */

    /**
     *  Calculate the orbital one-electron density matrix for a seniority-zero wave function expansion.
     * 
     *  @return The orbital (total, spin-summed) 1-DM.
     */
    template <typename Z = ONVBasis>
    enable_if_t<std::is_same<Z, SeniorityZeroONVBasis>::value, Orbital1DM<double>> calculate1DM() const {

        // Prepare some variables.
        const auto K = this->onv_basis.numberOfSpatialOrbitals();
        const auto dimension = this->onv_basis.dimension();
        Orbital1DM<double> D = Orbital1DM<double>::Zero(K);

        // Create the first ONV (with address 0). In DOCI, the ONV basis for alpha and beta is equal, so we can use the proxy ONV basis.
        const auto onv_basis_proxy = this->onv_basis.proxy();
        SpinUnresolvedONV onv = onv_basis_proxy.constructONVFromAddress(0);
        for (size_t I = 0; I < dimension; I++) {  // I loops over all the addresses of the doubly-occupied ONVs

            for (size_t e1 = 0; e1 < onv_basis_proxy.numberOfElectrons(); e1++) {  // e1 (electron 1) loops over the number of electrons
                const size_t p = onv.occupationIndexOf(e1);                        // retrieve the index of the orbital the electron occupies
                const double c_I = this->coefficient(I);                           // coefficient of the I-th basis vector

                D(p, p) += 2 * std::pow(c_I, 2);
            }

            if (I < dimension - 1) {  // prevent the last permutation from occurring
                onv_basis_proxy.transformONVToNextPermutation(onv);
            }
        }

        return D;
    }


    /**
     *  Calculate the two-electron density matrix for a seniority-zero wave function expansion.
     * 
     *  @return The orbital (total, spin-summed) 2-DM.
     */
    template <typename Z = ONVBasis>
    enable_if_t<std::is_same<Z, SeniorityZeroONVBasis>::value, Orbital2DM<double>> calculate2DM() const { return this->calculateSpinResolved2DM().orbitalDensity(); }

    /**
     *  Calculate the spin-resolved one-electron density matrix for a seniority-zero wave function expansion.
     * 
     *  @return The spin-resolved 1-DM.
     */
    template <typename Z = ONVBasis>
    enable_if_t<std::is_same<Z, SeniorityZeroONVBasis>::value, SpinResolved1DM<double>> calculateSpinResolved1DM() const { return SpinResolved1DM<double>::FromOrbital1DM(this->calculate1DM()); }


    /**
     *  Calculate the spin-resolved two-electron density matrix for a seniority-zero wave function expansion.
     * 
     *  @return The spin-resolved 2-DM.
     */
    template <typename Z = ONVBasis>
    enable_if_t<std::is_same<Z, SeniorityZeroONVBasis>::value, SpinResolved2DM<double>> calculateSpinResolved2DM() const {

        // Prepare some variables.
        const auto K = this->onv_basis.numberOfSpatialOrbitals();
        const auto dimension = this->onv_basis.dimension();


        // For seniority-zero linear expansions, we only have to calculate d_aaaa and d_aabb.
        PureSpinResolved2DMComponent<double> d_aaaa = PureSpinResolved2DMComponent<double>::Zero(K);
        MixedSpinResolved2DMComponent<double> d_aabb = MixedSpinResolved2DMComponent<double>::Zero(K);


        // Create the first ONV (with address 0). In DOCI, the ONV basis for alpha and beta is equal, so we can use the proxy ONV basis.
        const auto onv_basis_proxy = this->onv_basis.proxy();
        SpinUnresolvedONV onv = onv_basis_proxy.constructONVFromAddress(0);
        for (size_t I = 0; I < dimension; I++) {  // I loops over all the addresses of the spin strings
            for (size_t p = 0; p < K; p++) {      // p loops over SOs
                if (onv.annihilate(p)) {          // if p is occupied in I

                    const double c_I = this->coefficient(I);  // coefficient of the I-th basis vector
                    const double c_I_2 = std::pow(c_I, 2);    // square of c_I

                    d_aabb(p, p, p, p) += c_I_2;

                    for (size_t q = 0; q < p; q++) {                          // q loops over SOs with an index smaller than p
                        if (onv.create(q)) {                                  // if q is not occupied in I
                            const size_t J = onv_basis_proxy.addressOf(onv);  // the address of the coupling string
                            const double c_J = this->coefficient(J);          // coefficient of the J-th basis vector

                            d_aabb(p, q, p, q) += c_I * c_J;
                            d_aabb(q, p, q, p) += c_I * c_J;  // since we're looping for q < p

                            onv.annihilate(q);  // reset the spin string after previous creation on q
                        }

                        else {  // if q is occupied in I
                            d_aaaa(p, p, q, q) += c_I_2;
                            d_aaaa(q, q, p, p) += c_I_2;  // since we're looping for q < p

                            d_aaaa(p, q, q, p) -= c_I_2;
                            d_aaaa(q, p, p, q) -= c_I_2;  // since we're looping for q < p

                            d_aabb(p, p, q, q) += c_I_2;
                            d_aabb(q, q, p, p) += c_I_2;  // since we're looping for q < p
                        }
                    }
                    onv.create(p);  // reset the spin string after previous annihilation on p
                }
            }

            if (I < dimension - 1) {  // prevent the last permutation from occurring
                onv_basis_proxy.transformONVToNextPermutation(onv);
            }
        }

        // For seniority-zero linear expansions, we have additional symmetries (two_rdm_aaaa = two_rdm_bbbb, two_rdm_aabb = two_rdm_bbaa)
        return SpinResolved2DM<double> {d_aaaa, d_aabb, d_aabb, d_aaaa};
    }


    /*
     *  MARK: Density matrices for spin-resolved selected ONV bases
     */

    /**
     *  Calculate the spin-resolved one-electron density matrix for a spin-resolved selected wave function expansion.
     * 
     *  @return The spin-resolved 1-DM.
     */
    template <typename Z = ONVBasis>
    enable_if_t<std::is_same<Z, SpinResolvedSelectedONVBasis>::value, SpinResolved1DM<double>> calculateSpinResolved1DM() const {

        size_t K = this->onv_basis.numberOfOrbitals();
        size_t dim = onv_basis.dimension();

        SpinResolved1DMComponent<double> D_aa = SpinResolved1DMComponent<double>::Zero(K);
        SpinResolved1DMComponent<double> D_bb = SpinResolved1DMComponent<double>::Zero(K);


        for (size_t I = 0; I < dim; I++) {  // loop over all addresses (1)
            SpinResolvedONV configuration_I = this->onv_basis.onvWithIndex(I);
            SpinUnresolvedONV alpha_I = configuration_I.onv(Spin::alpha);
            SpinUnresolvedONV beta_I = configuration_I.onv(Spin::beta);

            double c_I = this->coefficient(I);


            // Calculate the diagonal of the 1-DMs
            for (size_t p = 0; p < K; p++) {

                if (alpha_I.isOccupied(p)) {
                    D_aa(p, p) += std::pow(c_I, 2);
                }

                if (beta_I.isOccupied(p)) {
                    D_bb(p, p) += std::pow(c_I, 2);
                }
            }


            // Calculate the off-diagonal elements, by going over all other ONVs
            for (size_t J = I + 1; J < dim; J++) {

                SpinResolvedONV configuration_J = this->onv_basis.onvWithIndex(J);
                SpinUnresolvedONV alpha_J = configuration_J.onv(Spin::alpha);
                SpinUnresolvedONV beta_J = configuration_J.onv(Spin::beta);

                double c_J = this->coefficient(J);


                // 1 electron excitation in alpha (i.e. 2 differences), 0 in beta
                if ((alpha_I.countNumberOfDifferences(alpha_J) == 2) && (beta_I.countNumberOfDifferences(beta_J) == 0)) {

                    // Find the orbitals that are occupied in one string, and aren't in the other
                    size_t p = alpha_I.findDifferentOccupations(alpha_J)[0];  // we're sure that there is only 1 element in the std::vector<size_t>
                    size_t q = alpha_J.findDifferentOccupations(alpha_I)[0];  // we're sure that there is only 1 element in the std::vector<size_t>

                    // Calculate the total sign, and include it in the DM contribution
                    int sign = alpha_I.operatorPhaseFactor(p) * alpha_J.operatorPhaseFactor(q);
                    D_aa(p, q) += sign * c_I * c_J;
                    D_aa(q, p) += sign * c_I * c_J;
                }


                // 1 electron excitation in beta, 0 in alpha
                if ((alpha_I.countNumberOfDifferences(alpha_J) == 0) && (beta_I.countNumberOfDifferences(beta_J) == 2)) {

                    // Find the orbitals that are occupied in one string, and aren't in the other
                    size_t p = beta_I.findDifferentOccupations(beta_J)[0];  // we're sure that there is only 1 element in the std::vector<size_t>
                    size_t q = beta_J.findDifferentOccupations(beta_I)[0];  // we're sure that there is only 1 element in the std::vector<size_t>

                    // Calculate the total sign, and include it in the DM contribution
                    int sign = beta_I.operatorPhaseFactor(p) * beta_J.operatorPhaseFactor(q);
                    D_bb(p, q) += sign * c_I * c_J;
                    D_bb(q, p) += sign * c_I * c_J;
                }

            }  // loop over addresses J > I
        }      // loop over addresses I

        return SpinResolved1DM<double>(D_aa, D_bb);  // the total 1-DM is the sum of the spin components
    }


    /**
     *  Calculate the spin-resolved two-electron density matrix for a spin-resolved selected wave function expansion.
     * 
     *  @return The spin-resolved 2-DM.
     */
    template <typename Z = ONVBasis>
    enable_if_t<std::is_same<Z, SpinResolvedSelectedONVBasis>::value, SpinResolved2DM<double>> calculateSpinResolved2DM() const {

        size_t K = this->onv_basis.numberOfOrbitals();
        size_t dim = onv_basis.dimension();

        PureSpinResolved2DMComponent<double> d_aaaa = PureSpinResolved2DMComponent<double>::Zero(K);
        MixedSpinResolved2DMComponent<double> d_aabb = MixedSpinResolved2DMComponent<double>::Zero(K);
        MixedSpinResolved2DMComponent<double> d_bbaa = MixedSpinResolved2DMComponent<double>::Zero(K);
        PureSpinResolved2DMComponent<double> d_bbbb = PureSpinResolved2DMComponent<double>::Zero(K);

        for (size_t I = 0; I < dim; I++) {  // loop over all addresses I

            SpinResolvedONV configuration_I = this->onv_basis.onvWithIndex(I);
            SpinUnresolvedONV alpha_I = configuration_I.onv(Spin::alpha);
            SpinUnresolvedONV beta_I = configuration_I.onv(Spin::beta);

            double c_I = this->coefficient(I);

            for (size_t p = 0; p < K; p++) {

                // 'Diagonal' elements of the 2-DM: aaaa and aabb
                if (alpha_I.isOccupied(p)) {
                    for (size_t q = 0; q < K; q++) {
                        if (beta_I.isOccupied(q)) {
                            d_aabb(p, p, q, q) += std::pow(c_I, 2);
                        }

                        if (p != q) {  // can't create/annihilate the same orbital twice
                            if (alpha_I.isOccupied(q)) {
                                d_aaaa(p, p, q, q) += std::pow(c_I, 2);
                                d_aaaa(p, q, q, p) -= std::pow(c_I, 2);
                            }
                        }

                    }  // loop over q
                }

                // 'Diagonal' elements of the 2-DM: bbbb and bbaa
                if (beta_I.isOccupied(p)) {
                    for (size_t q = 0; q < K; q++) {
                        if (alpha_I.isOccupied(q)) {
                            d_bbaa(p, p, q, q) += std::pow(c_I, 2);
                        }

                        if (p != q) {  // can't create/annihilate the same orbital twice
                            if (beta_I.isOccupied(q)) {
                                d_bbbb(p, p, q, q) += std::pow(c_I, 2);
                                d_bbbb(p, q, q, p) -= std::pow(c_I, 2);
                            }
                        }
                    }  // loop over q
                }
            }  // loop over q


            for (size_t J = I + 1; J < dim; J++) {

                SpinResolvedONV configuration_J = this->onv_basis.onvWithIndex(J);
                SpinUnresolvedONV alpha_J = configuration_J.onv(Spin::alpha);
                SpinUnresolvedONV beta_J = configuration_J.onv(Spin::beta);

                double c_J = this->coefficient(J);

                // 1 electron excitation in alpha, 0 in beta
                if ((alpha_I.countNumberOfDifferences(alpha_J) == 2) && (beta_I.countNumberOfDifferences(beta_J) == 0)) {

                    // Find the orbitals that are occupied in one string, and aren't in the other
                    size_t p = alpha_I.findDifferentOccupations(alpha_J)[0];  // we're sure that there is only 1 element in the std::vector<size_t>
                    size_t q = alpha_J.findDifferentOccupations(alpha_I)[0];  // we're sure that there is only 1 element in the std::vector<size_t>

                    // Calculate the total sign
                    int sign = alpha_I.operatorPhaseFactor(p) * alpha_J.operatorPhaseFactor(q);


                    for (size_t r = 0; r < K; r++) {  // r loops over spatial orbitals

                        if (alpha_I.isOccupied(r) && alpha_J.isOccupied(r)) {  // r must be occupied on the left and on the right
                            if ((p != r) && (q != r)) {                        // can't create or annihilate the same orbital
                                // Fill in the 2-DM contributions
                                d_aaaa(p, q, r, r) += sign * c_I * c_J;
                                d_aaaa(r, q, p, r) -= sign * c_I * c_J;
                                d_aaaa(p, r, r, q) -= sign * c_I * c_J;
                                d_aaaa(r, r, p, q) += sign * c_I * c_J;

                                d_aaaa(q, p, r, r) += sign * c_I * c_J;
                                d_aaaa(q, r, r, p) -= sign * c_I * c_J;
                                d_aaaa(r, p, q, r) -= sign * c_I * c_J;
                                d_aaaa(r, r, q, p) += sign * c_I * c_J;
                            }
                        }

                        if (beta_I.isOccupied(r)) {  // beta_I == beta_J from the previous if-branch

                            // Fill in the 2-DM contributions
                            d_aabb(p, q, r, r) += sign * c_I * c_J;
                            d_aabb(q, p, r, r) += sign * c_I * c_J;

                            d_bbaa(r, r, p, q) += sign * c_I * c_J;
                            d_bbaa(r, r, q, p) += sign * c_I * c_J;
                        }
                    }
                }


                // 0 electron excitations in alpha, 1 in beta
                if ((alpha_I.countNumberOfDifferences(alpha_J) == 0) && (beta_I.countNumberOfDifferences(beta_J) == 2)) {

                    // Find the orbitals that are occupied in one string, and aren't in the other
                    size_t p = beta_I.findDifferentOccupations(beta_J)[0];  // we're sure that there is only 1 element in the std::vector<size_t>
                    size_t q = beta_J.findDifferentOccupations(beta_I)[0];  // we're sure that there is only 1 element in the std::vector<size_t>

                    // Calculate the total sign
                    int sign = beta_I.operatorPhaseFactor(p) * beta_J.operatorPhaseFactor(q);


                    for (size_t r = 0; r < K; r++) {  // r loops over spatial orbitals

                        if (beta_I.isOccupied(r) && beta_J.isOccupied(r)) {  // r must be occupied on the left and on the right
                            if ((p != r) && (q != r)) {                      // can't create or annihilate the same orbital
                                // Fill in the 2-DM contributions
                                d_bbbb(p, q, r, r) += sign * c_I * c_J;
                                d_bbbb(r, q, p, r) -= sign * c_I * c_J;
                                d_bbbb(p, r, r, q) -= sign * c_I * c_J;
                                d_bbbb(r, r, p, q) += sign * c_I * c_J;

                                d_bbbb(q, p, r, r) += sign * c_I * c_J;
                                d_bbbb(q, r, r, p) -= sign * c_I * c_J;
                                d_bbbb(r, p, q, r) -= sign * c_I * c_J;
                                d_bbbb(r, r, q, p) += sign * c_I * c_J;
                            }
                        }

                        if (alpha_I.isOccupied(r)) {  // alpha_I == alpha_J from the previous if-branch

                            // Fill in the 2-DM contributions
                            d_bbaa(p, q, r, r) += sign * c_I * c_J;
                            d_bbaa(q, p, r, r) += sign * c_I * c_J;

                            d_aabb(r, r, p, q) += sign * c_I * c_J;
                            d_aabb(r, r, q, p) += sign * c_I * c_J;
                        }
                    }
                }


                // 1 electron excitation in alpha, 1 in beta
                if ((alpha_I.countNumberOfDifferences(alpha_J) == 2) && (beta_I.countNumberOfDifferences(beta_J) == 2)) {

                    // Find the orbitals that are occupied in one string, and aren't in the other
                    size_t p = alpha_I.findDifferentOccupations(alpha_J)[0];  // we're sure that there is only 1 element in the std::vector<size_t>
                    size_t q = alpha_J.findDifferentOccupations(alpha_I)[0];  // we're sure that there is only 1 element in the std::vector<size_t>

                    size_t r = beta_I.findDifferentOccupations(beta_J)[0];  // we're sure that there is only 1 element in the std::vector<size_t>
                    size_t s = beta_J.findDifferentOccupations(beta_I)[0];  // we're sure that there is only 1 element in the std::vector<size_t>

                    // Calculate the total sign, and include it in the 2-DM contribution
                    int sign = alpha_I.operatorPhaseFactor(p) * alpha_J.operatorPhaseFactor(q) * beta_I.operatorPhaseFactor(r) * beta_J.operatorPhaseFactor(s);
                    d_aabb(p, q, r, s) += sign * c_I * c_J;
                    d_aabb(q, p, s, r) += sign * c_I * c_J;

                    d_bbaa(r, s, p, q) += sign * c_I * c_J;
                    d_bbaa(s, r, q, p) += sign * c_I * c_J;
                }


                // 2 electron excitations in alpha, 0 in beta
                if ((alpha_I.countNumberOfDifferences(alpha_J) == 4) && (beta_I.countNumberOfDifferences(beta_J) == 0)) {

                    // Find the orbitals that are occupied in one string, and aren't in the other
                    std::vector<size_t> occupied_indices_I = alpha_I.findDifferentOccupations(alpha_J);  // we're sure this has two elements
                    size_t p = occupied_indices_I[0];
                    size_t r = occupied_indices_I[1];

                    std::vector<size_t> occupied_indices_J = alpha_J.findDifferentOccupations(alpha_I);  // we're sure this has two elements
                    size_t q = occupied_indices_J[0];
                    size_t s = occupied_indices_J[1];


                    // Calculate the total sign, and include it in the 2-DM contribution
                    int sign = alpha_I.operatorPhaseFactor(p) * alpha_I.operatorPhaseFactor(r) * alpha_J.operatorPhaseFactor(q) * alpha_J.operatorPhaseFactor(s);
                    d_aaaa(p, q, r, s) += sign * c_I * c_J;
                    d_aaaa(p, s, r, q) -= sign * c_I * c_J;
                    d_aaaa(r, q, p, s) -= sign * c_I * c_J;
                    d_aaaa(r, s, p, q) += sign * c_I * c_J;

                    d_aaaa(q, p, s, r) += sign * c_I * c_J;
                    d_aaaa(s, p, q, r) -= sign * c_I * c_J;
                    d_aaaa(q, r, s, p) -= sign * c_I * c_J;
                    d_aaaa(s, r, q, p) += sign * c_I * c_J;
                }


                // 0 electron excitations in alpha, 2 in beta
                if ((alpha_I.countNumberOfDifferences(alpha_J) == 0) && (beta_I.countNumberOfDifferences(beta_J) == 4)) {

                    // Find the orbitals that are occupied in one string, and aren't in the other
                    std::vector<size_t> occupied_indices_I = beta_I.findDifferentOccupations(beta_J);  // we're sure this has two elements
                    size_t p = occupied_indices_I[0];
                    size_t r = occupied_indices_I[1];

                    std::vector<size_t> occupied_indices_J = beta_J.findDifferentOccupations(beta_I);  // we're sure this has two elements
                    size_t q = occupied_indices_J[0];
                    size_t s = occupied_indices_J[1];


                    // Calculate the total sign, and include it in the 2-DM contribution
                    int sign = beta_I.operatorPhaseFactor(p) * beta_I.operatorPhaseFactor(r) * beta_J.operatorPhaseFactor(q) * beta_J.operatorPhaseFactor(s);
                    d_bbbb(p, q, r, s) += sign * c_I * c_J;
                    d_bbbb(p, s, r, q) -= sign * c_I * c_J;
                    d_bbbb(r, q, p, s) -= sign * c_I * c_J;
                    d_bbbb(r, s, p, q) += sign * c_I * c_J;

                    d_bbbb(q, p, s, r) += sign * c_I * c_J;
                    d_bbbb(s, p, q, r) -= sign * c_I * c_J;
                    d_bbbb(q, r, s, p) -= sign * c_I * c_J;
                    d_bbbb(s, r, q, p) += sign * c_I * c_J;
                }

            }  // loop over all addresses J > I

        }  // loop over all addresses I

        return SpinResolved2DM<double>(d_aaaa, d_aabb, d_bbaa, d_bbbb);
    }


    /**
     *  Calculate the one-electron density matrix for a spin-resolved selected wave function expansion.
     * 
     *  @return The orbital (total, spin-summed) 1-DM.
     */
    template <typename Z = ONVBasis>
    enable_if_t<std::is_same<Z, SpinResolvedSelectedONVBasis>::value, Orbital1DM<double>> calculate1DM() const { return this->calculateSpinResolved1DM().orbitalDensity(); }


    /**
     *  Calculate the two-electron density matrix for a spin-resolved selected wave function expansion.
     * 
     *  @return The orbital (total, spin-summed) 2-DM.
     */
    template <typename Z = ONVBasis>
    enable_if_t<std::is_same<Z, SpinResolvedSelectedONVBasis>::value, Orbital2DM<double>> calculate2DM() const { return this->calculateSpinResolved2DM().orbitalDensity(); }


    /**
     *  MARK: Entropy
     */

    /**
     *  @return The Shannon entropy (information content) of the wave function.
     */
    double calculateShannonEntropy() const {

        // Sum over the ONV basis dimension, and only include the term if c_k != 0
        // We might as well replace all coefficients that are 0 by 1, since log(1) = 0 so there is no influence on the final entropy value
        Eigen::ArrayXd coefficients_replaced = this->m_coefficients.unaryExpr([](double c) { return c < 1.0e-18 ? 1 : c; });  // replace 0 by 1

        Eigen::ArrayXd coefficients_squared = coefficients_replaced.square();
        Eigen::ArrayXd log_coefficients_squared = coefficients_squared.log();  // natural logarithm (ln)

        return -1 / std::log(2) * (coefficients_squared * log_coefficients_squared).sum();
    }


    /*
     *  MARK: Iterating
     */

    /**
     *  Iterate over all expansion coefficients and corresponding ONVs, and apply the given callback function.
     * 
     *  @param callback                 The function to be applied in every iteration. Its arguments are an expansion coefficient and the corresponding ONV.
     */
    template <typename Z = ONVBasis>
    enable_if_t<std::is_same<Z, SpinResolvedONVBasis>::value, void> forEach(const std::function<void(const double, const SpinResolvedONV)>& callback) const {

        // Iterate over all ONVs in this ONV basis, and look up the corresponding coefficient.
        const auto& onv_basis = this->onv_basis;
        const auto& coefficients = this->m_coefficients;
        onv_basis.forEach([&onv_basis, &callback, &coefficients](const SpinUnresolvedONV& onv_alpha, const size_t I_alpha, const SpinUnresolvedONV& onv_beta, const size_t I_beta) {
            const SpinResolvedONV onv {onv_alpha, onv_beta};
            const auto address = onv_basis.compoundAddress(I_alpha, I_beta);
            const auto coefficient = coefficients(address);

            callback(coefficient, onv);
        });
    }


    /*
     *  MARK: Comparing
     */

    /** 
     *  @param other            wave function for the comparison
     *  @param tolerance        tolerance for the comparison of coefficients
     * 
     *  @return if two wave functions are equal within a given tolerance
     */
    bool isApprox(const LinearExpansion<ONVBasis>& other, double tolerance = 1e-10) const {

        if (this->onv_basis.dimension() != other.onv_basis.dimension()) {
            return false;
        }

        return (this->coefficients()).isEqualEigenvectorAs(other.coefficients(), tolerance);
    }
};


}  // namespace GQCP
