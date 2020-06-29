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


#include "Basis/Integrals/BaseOneElectronIntegralEngine.hpp"
#include "Basis/Integrals/BaseTwoElectronIntegralEngine.hpp"
#include "Basis/Integrals/IntegralEngine.hpp"
#include "Basis/Integrals/Interfaces/LibcintInterfacer.hpp"
#include "Basis/Integrals/Interfaces/LibintInterfacer.hpp"
#include "Basis/ScalarBasis/ScalarBasis.hpp"
#include "Basis/ScalarBasis/ShellSet.hpp"
#include "Mathematical/Representation/QCMatrix.hpp"
#include "Mathematical/Representation/QCRankFourTensor.hpp"
#include "Operator/FirstQuantized/Operator.hpp"
#include "Utilities/aliases.hpp"

#include <algorithm>
#include <array>
#include <memory>


namespace GQCP {


/**
 *  A class that calculates integrals over ShellSets: it loops over all shells in the given ShellSets. It just has two static member functions, for cleaner interface calls.
 */
class IntegralCalculator {
public:
    /*
     *  PUBLIC METHODS
     */

    /**
     *  Calculate all one-electron integrals over the basis functions inside the given shell sets.
     * 
     *  @param engine                   the engine that can calculate one-electron integrals over shells (not const because we allow for non-const Engine::calculate() calls)
     *  @param left_shell_set           the set of shells that should appear on the left of the operator
     *  @param right_shell_set          the set of shells that should appear on the right of the operator
     * 
     *  @tparam Shell                   the type of shell the integral engine is able to handle
     *  @tparam N                       the number of components the operator has
     *  @tparam IntegralScalar          the scalar representation of an integral
     */
    template <typename Shell, size_t N, typename IntegralScalar>
    static auto calculate(BaseOneElectronIntegralEngine<Shell, N, IntegralScalar>& engine, const ShellSet<Shell>& left_shell_set, const ShellSet<Shell>& right_shell_set) -> std::array<MatrixX<IntegralScalar>, N> {

        // Initialize the N components of the matrix representations of the operator.
        const auto nbf_left = left_shell_set.numberOfBasisFunctions();
        const auto nbf_right = right_shell_set.numberOfBasisFunctions();

        std::array<MatrixX<IntegralScalar>, N> components;
        for (auto& component : components) {
            component = MatrixX<IntegralScalar>::Zero(nbf_left, nbf_right);
        }


        // Loop over all left and right shells and let the engine calculate the integrals over the pairs of shells.
        const auto nsh_left = left_shell_set.numberOfShells();
        const auto left_shells = left_shell_set.asVector();
        const auto nsh_right = right_shell_set.numberOfShells();
        const auto right_shells = right_shell_set.asVector();

        for (size_t left_shell_index = 0; left_shell_index < nsh_left; left_shell_index++) {
            const auto left_bf_index = left_shell_set.basisFunctionIndex(left_shell_index);
            const auto left_shell = left_shells[left_shell_index];

            for (size_t right_shell_index = 0; right_shell_index < nsh_right; right_shell_index++) {
                const auto right_bf_index = right_shell_set.basisFunctionIndex(right_shell_index);
                const auto right_shell = right_shells[right_shell_index];

                const auto buffer = engine.calculate(left_shell, right_shell);

                // Only if the integrals are not all zero, place them inside the full matrices.
                if (buffer->areIntegralsAllZero()) {
                    continue;
                }
                buffer->emplace(components, left_bf_index, right_bf_index);
            }  // right shells loop
        }      // left shells loop

        return components;
    }


    /**
     *  Calculate all two-electron integrals over the basis functions inside the given shell sets.
     * 
     *  @param engine                       the engine that can calculate two-electron integrals over shells
     *  @param left_shell_set               the set of shells that should appear on the left of the operator
     *  @param right_shell_set              the set of shells that should appear on the right of the operator
     * 
     *  @tparam Shell                       the type of shell the integral engine is able to handle
     *  @tparam N                           the number of components the operator has
     *  @tparam IntegralScalar              the scalar representation of an integral
     */
    template <typename Shell, size_t N, typename IntegralScalar>
    static auto calculate(BaseTwoElectronIntegralEngine<Shell, N, IntegralScalar>& engine, const ShellSet<Shell>& left_shell_set, const ShellSet<Shell>& right_shell_set) -> std::array<Tensor<IntegralScalar, 4>, N> {

        // Initialize the N components of the matrix representations of the operator.
        const auto nbf_left = left_shell_set.numberOfBasisFunctions();
        const auto nbf_right = right_shell_set.numberOfBasisFunctions();

        std::array<Tensor<IntegralScalar, 4>, N> components;
        for (auto& component : components) {
            component = Tensor<IntegralScalar, 4>(nbf_left, nbf_left, nbf_right, nbf_right);
            component.setZero();
        }


        // Loop over all left and right shells and let the engine calculate the integrals over the 4-tuple of shells.
        const auto nsh_left = left_shell_set.numberOfShells();
        const auto left_shells = left_shell_set.asVector();
        const auto nsh_right = right_shell_set.numberOfShells();
        const auto right_shells = right_shell_set.asVector();

        for (size_t left_shell_index1 = 0; left_shell_index1 < nsh_left; left_shell_index1++) {
            const auto left_bf1_index = left_shell_set.basisFunctionIndex(left_shell_index1);
            const auto left_shell1 = left_shells[left_shell_index1];

            for (size_t left_shell_index2 = 0; left_shell_index2 < nsh_left; left_shell_index2++) {
                const auto left_bf2_index = left_shell_set.basisFunctionIndex(left_shell_index2);
                const auto left_shell2 = left_shells[left_shell_index2];

                for (size_t right_shell_index1 = 0; right_shell_index1 < nsh_right; right_shell_index1++) {
                    const auto right_bf1_index = right_shell_set.basisFunctionIndex(right_shell_index1);
                    const auto right_shell1 = right_shells[right_shell_index1];

                    for (size_t right_shell_index2 = 0; right_shell_index2 < nsh_right; right_shell_index2++) {
                        const auto right_bf2_index = right_shell_set.basisFunctionIndex(right_shell_index2);
                        const auto right_shell2 = right_shells[right_shell_index2];

                        const auto buffer = engine.calculate(left_shell1, left_shell2, right_shell1, right_shell2);

                        // Only if the integrals are not all zero, place them inside the full matrices
                        if (buffer->areIntegralsAllZero()) {
                            continue;
                        }
                        buffer->emplace(components, left_bf1_index, left_bf2_index, right_bf1_index, right_bf2_index);  // place the calculated integrals inside the full tensors

                    }  // sh4_index
                }      // right_shell_index1
            }          // left_shell_index2
        }              // left_shell_index1

        return components;
    }


    /*
     *  PUBLIC METHODS - LIBINT2 INTEGRALS
     */

    /**
     *  Calculate the integrals over the given (first-quantized) one-electron operator, over a left and right scalar basis, using Libint2.
     * 
     *  @param fq_one_op                            the first-quantized one-electron operator
     *  @param left_scalar_basis                    the scalar basis that contains the shells that appear to the left of the operator
     *  @param right_scalar_basis                   the scalar basis that contains the shell sets that appear to the right of the operator over which the integrals should be calculated
     * 
     *  @tparam FQOneElectronOperator               the type of the first-quantized one-electron operator
     * 
     *  @return the matrix representation (integrals) of the given first-quantized operator in this scalar basis
     */
    template <typename FQOneElectronOperator>
    static Matrix<double> calculateLibintIntegrals(const FQOneElectronOperator& fq_one_op, const ScalarBasis<GTOShell>& left_scalar_basis, const ScalarBasis<GTOShell>& right_scalar_basis) {

        const auto left_shell_set = left_scalar_basis.shellSet();
        const auto right_shell_set = right_scalar_basis.shellSet();

        // Construct the libint engine
        const auto max_nprim = std::max(left_shell_set.maximumNumberOfPrimitives(), right_shell_set.maximumNumberOfPrimitives());
        const auto max_l = std::max(left_shell_set.maximumAngularMomentum(), right_shell_set.maximumAngularMomentum());
        auto engine = IntegralEngine::Libint(fq_one_op, max_nprim, max_l);


        // Calculate the integrals using the engine
        const auto integrals = IntegralCalculator::calculate(engine, left_shell_set, right_shell_set);
        return integrals[0];
    }


    /**
     *  Calculate the integrals over the given (first-quantized) one-electron operator, within a given scalar basis, using Libint2.
     * 
     *  @param fq_one_op                            the first-quantized one-electron operator
     *  @param scalar_basis                         the scalar basis that contains the shells over which the integrals should be calculated
     * 
     *  @tparam FQOneElectronOperator               the type of the first-quantized one-electron operator
     * 
     *  @return the matrix representation (integrals) of the given first-quantized operator in this scalar basis
     */
    template <typename FQOneElectronOperator>
    static QCMatrix<double> calculateLibintIntegrals(const FQOneElectronOperator& fq_one_op, const ScalarBasis<GTOShell>& scalar_basis) {

        return QCMatrix<double>(IntegralCalculator::calculateLibintIntegrals(fq_one_op, scalar_basis, scalar_basis));  // the same scalar basis appear on the left and right of the operator
    }


    /**
     *  Calculate the integrals over the given (first-quantized) one-electron operator, over a left and right scalar basis, using Libint2.
     * 
     *  @param fq_one_op                            the first-quantized one-electron operator
     *  @param left_scalar_basis                    the scalar basis that contains the shells that should appear to the left of the operator
     *  @param right_scalar_basis                   the scalar basis that contains the shells that should appear to the right of the operator
     * 
     *  @return the matrix representation (integrals) of the given first-quantized operator in this scalar basis
     */
    static std::array<Matrix<double>, 3> calculateLibintIntegrals(const ElectronicDipoleOperator& fq_one_op, const ScalarBasis<GTOShell>& left_scalar_basis, const ScalarBasis<GTOShell>& right_scalar_basis) {

        const auto left_shell_set = left_scalar_basis.shellSet();
        const auto right_shell_set = right_scalar_basis.shellSet();

        // Construct the libint engine
        const auto max_nprim = std::max(left_shell_set.maximumNumberOfPrimitives(), right_shell_set.maximumNumberOfPrimitives());
        const auto max_l = std::max(left_shell_set.maximumAngularMomentum(), right_shell_set.maximumAngularMomentum());
        auto engine = IntegralEngine::Libint(fq_one_op, max_nprim, max_l);


        // Calculate the integrals using the engine
        const auto integrals = IntegralCalculator::calculate(engine, left_shell_set, right_shell_set);
        return {integrals[0], integrals[1], integrals[2]};
    }


    /**
     *  Calculate the integrals over the given (first-quantized) one-electron operator, within a given scalar basis, using Libint2.
     * 
     *  @param fq_one_op                            the first-quantized one-electron operator
     *  @param scalar_basis                         the scalar basis that contains the shells over which the integrals should be calculated
     * 
     *  @return the matrix representation (integrals) of the given first-quantized operator in this scalar basis
     */
    static std::array<QCMatrix<double>, 3> calculateLibintIntegrals(const ElectronicDipoleOperator& fq_one_op, const ScalarBasis<GTOShell>& scalar_basis) {

        // Convert the array of Matrix to an array of QCMatrix
        const auto calculated_components = IntegralCalculator::calculateLibintIntegrals(fq_one_op, scalar_basis, scalar_basis);  // the same scalar basis appear on the left and right of the operator

        std::array<QCMatrix<double>, 3> converted_components {};
        for (size_t i = 0; i < 3; i++) {
            converted_components[i] = QCMatrix<double>(calculated_components[i]);
        }
        return converted_components;
    }


    /**
     *  Calculate the integrals over the given (first-quantized) two-electron operator, within a given scalar basis, using Libint2.
     * 
     *  @param fq_two_op                    the first-quantized operator
     *  @param scalar_basis                 the scalar basis that contains the shells over which the integrals should be calculated
     * 
     *  @return the matrix representation (integrals) of the given first-quantized operator in this scalar basis
     */
    static QCRankFourTensor<double> calculateLibintIntegrals(const CoulombRepulsionOperator& fq_two_op, const ScalarBasis<GTOShell>& scalar_basis) {

        return QCRankFourTensor<double>(IntegralCalculator::calculateLibintIntegrals(fq_two_op, scalar_basis, scalar_basis));  // the same scalar basis appear on the left and right of the operator
    }


    /**
     *  Calculate the integrals over the given (first-quantized) two-electron operator, over a left and right scalar basis, using Libint2.
     * 
     *  @param fq_two_op                            the first-quantized operator
     *  @param left_scalar_basis                    the scalar basis that contains the shells that should appear to the left of the operator
     *  @param right_scalar_basis                   the scalar basis that contains the shells that should appear to the right of the operator
     * 
     *  @return the matrix representation (integrals) of the given first-quantized operator in this scalar basis
     */
    static Tensor<double, 4> calculateLibintIntegrals(const CoulombRepulsionOperator& fq_two_op, const ScalarBasis<GTOShell>& left_scalar_basis, const ScalarBasis<GTOShell>& right_scalar_basis) {

        const auto left_shell_set = left_scalar_basis.shellSet();
        const auto right_shell_set = right_scalar_basis.shellSet();

        // Construct the libint engine
        const auto max_nprim = std::max(left_shell_set.maximumNumberOfPrimitives(), right_shell_set.maximumNumberOfPrimitives());
        const auto max_l = std::max(left_shell_set.maximumAngularMomentum(), right_shell_set.maximumAngularMomentum());
        auto engine = IntegralEngine::Libint(fq_two_op, max_nprim, max_l);


        // Calculate the integrals using the engine
        const auto integrals = IntegralCalculator::calculate(engine, left_shell_set, right_shell_set);
        return integrals[0];
    }


    /*
     *  PUBLIC METHODS - LIBCINT INTEGRALS
     *  Note that the Libcint integrals should only be used for Cartesian ShellSets
     */

    /**
     *  Calculate the integrals over the given (first-quantized) one-electron operator, within a given scalar basis, using libcint.
     *
     *  @param fq_one_op                            the first-quantized one-electron operator
     *  @param scalar_basis                         the scalar basis that contains the shells over which the integrals should be calculated
     * 
     *  @tparam FQOneElectronOperator               the type of the first-quantized one-electron operator
     * 
     *  @note Only use this function for all-Cartesian ShellSets.
     * 
     *  @return the matrix representation of the overlap operator in this AO basis, using the libcint integral engine
     */
    template <typename FQOneElectronOperator>
    static QCMatrix<double> calculateLibcintIntegrals(const FQOneElectronOperator& fq_one_op, const ScalarBasis<GTOShell>& scalar_basis) {

        const auto shell_set = scalar_basis.shellSet();

        auto engine = IntegralEngine::Libcint(fq_one_op, shell_set);
        const auto integrals = IntegralCalculator::calculate(engine, shell_set, shell_set);
        return integrals[0];
    }


    /**
     *  Calculate the integrals over the given electronic dipole operator, within a given scalar basis, using libcint.
     *
     *  @param fq_one_op                            the electronic dipole operator
     *  @param scalar_basis                         the scalar basis that contains the shells over which the integrals should be calculated
     *
     *  @note Only use this function for all-Cartesian ShellSets.
     * 
     *  @return the matrix representation of the Cartesian components of the electrical dipole operator in this AO basis, using the libcint integral engine
     */
    static std::array<QCMatrix<double>, 3> calculateLibcintIntegrals(const ElectronicDipoleOperator& fq_one_op, const ScalarBasis<GTOShell>& scalar_basis) {

        const auto shell_set = scalar_basis.shellSet();

        auto engine = IntegralEngine::Libcint(fq_one_op, shell_set);
        const auto integrals = IntegralCalculator::calculate(engine, shell_set, shell_set);
        return {integrals[0], integrals[1], integrals[2]};
    }


    /**
     *  Calculate the Coulomb repulsion energy integrals, within a given scalar basis, using Libcint.
     *
     *  @param fq_op                                the first-quantized operator
     *  @param scalar_basis                         the scalar basis that contains the shells over which the integrals should be calculated
     * 
     *  @note Only use this function for all-Cartesian ShellSets.
     * 
     *  @return the matrix representation of the Coulomb repulsion operator in this AO basis, using the libcint integral engine
     */
    static QCRankFourTensor<double> calculateLibcintIntegrals(const CoulombRepulsionOperator& fq_op, const ScalarBasis<GTOShell>& scalar_basis) {

        const auto shell_set = scalar_basis.shellSet();

        auto engine = IntegralEngine::Libcint(fq_op, shell_set);
        const auto integrals = IntegralCalculator::calculate(engine, shell_set, shell_set);
        return integrals[0];
    }
};


}  // namespace GQCP
