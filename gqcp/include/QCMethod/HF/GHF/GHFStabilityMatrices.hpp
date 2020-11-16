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


#include "Mathematical/Representation/ImplicitRankFourTensorSlice.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/QCStructure.hpp"
#include "Utilities/aliases.hpp"


namespace GQCP {
namespace QCMethod {


/**
 *  The generalized Hartree-Fock quantum chemical method.
 * 
 *  @tparam _Scalar             the type of scalar that is used for the expansion of the spinors in their underlying scalar basis
 */
template <typename _Scalar>
class GHFStability {

public:
    using Scalar = _Scalar;


public:
    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return the partial stability matrix `A` from the GHF stability conditions.
     * 
     *  @note The formula for the `A` matrix is as follows:
     *      A_IAJB = \delta_IJ * F_BA - \delta_BA * F_IJ + (AI||JB)
     * 
     *  @param ghf_structure        The GHF QCStructure which contains the ground state parameters of the performed GHF calculation.
     *  @param gsq_hamiltonian      The generalised, second quantized hamiltonian, which contains the necessary two electron operators.
     */
    const GQCP::Matrix<Scalar> calculatePartialStabilityMatrixA(const QCStructure<QCModel::GHF<Scalar>>& ghf_structure, const GSQHamiltonian<Scalar>& gsq_hamiltonian) const {

        // Get the ground state parameters from the given QCStructure.
        const auto& parameters = ghf_structure.groundStateParameters();

        // Determine the number of occupied and virtual orbitals.
        const auto& number_of_occupied_orbitals = parameters.numberOfElectrons();
        const auto& number_of_virtual_orbitals = gsq_hamiltonian.numberOfOrbitals() - number_of_occupied_orbitals;

        // We need the two-electron integrals in MO basis, hence why we transform them with the coefficient matrix.
        // The ground state coefficient matrix is obtained from the QCModel.
        // We need the anti-symmetrized tensor: (AI||JB) = (AI|JB) - (AB|JI). This is obtained by the `.antisymmetrized()` method.
        const auto& g = gsq_hamiltonian.twoElectron().transformed(parameters.coefficientMatrix()).antisymmetrized();

        // The elements F_BA and F_IJ are the eigenvalues of the one-electron Fock operator. This can be done using the MO energies.
        // The MO energies are contained within the QC Structure.
        std::vector<Scalar> mo_energies;  // We use a std::vector in order to be able to slice the vector later on.

        for (int i = 0; i < parameters.numberOfSpinors(); i++) {
            mo_energies.push_back(parameters.orbitalEnergy(i));
        }

        // The MO energies should be split in the respective real and virtual parts.
        std::vector<Scalar> mo_energies_occupied;
        std::copy(mo_energies.begin(), mo_energies.begin() + number_of_occupied_orbitals, std::back_inserter(mo_energies_occupied));

        std::vector<Scalar> mo_energies_virtual;
        std::copy(mo_energies.begin() + number_of_occupied_orbitals, mo_energies.end(), std::back_inserter(mo_energies_virtual));

        // We create a matrix containing the correct elements needed in the formula.
        // These values will be added to the final tensor in the last step of this calculation.
        GQCP::MatrixX<Scalar> F_values(number_of_virtual_orbitals, number_of_occupied_orbitals);

        for (int a = 0; a < number_of_virtual_orbitals; a++) {
            for (int i = 0; i < number_of_occupied_orbitals; i++) {
                F_values(a, i) = mo_energies_virtual[a] + (-1 * mo_energies_occupied[i]);
            }
        }

        // The next step is to create the needed tensor slice.
        GQCP::Tensor<Scalar, 4> A_iajb(number_of_occupied_orbitals, number_of_virtual_orbitals, number_of_occupied_orbitals, number_of_virtual_orbitals);
        for (int i = 0; i < number_of_occupied_orbitals; i++) {
            for (int a = number_of_occupied_orbitals; a < parameters.numberOfSpinors(); a++) {
                for (int j = 0; j < number_of_occupied_orbitals; j++) {
                    for (int b = number_of_occupied_orbitals; b < parameters.numberOfSpinors(); b++) {
                        A_iajb(i, a - number_of_occupied_orbitals, j, b - number_of_occupied_orbitals) = g.parameters()(a, i, j, b);
                    }
                }
            }
        }

        // Add the previously calculated F elements on the correct positions.
        for (int a = 0; a < number_of_virtual_orbitals; a++) {
            for (int i = 0; i < number_of_occupied_orbitals; i++) {
                A_iajb(i, a, i, a) += F_values(a, i);
            }
        }
        A_iajb.print();
        //Finally, reshape the tensor to a matrix.
        const auto A_matrix = A_iajb.reshape(number_of_occupied_orbitals * number_of_virtual_orbitals, number_of_occupied_orbitals * number_of_virtual_orbitals);

        return A_matrix;
    }


    /**
     *  @return the partial stability matrix `B` from the GHF stability conditions.
     *
     *  @note The formula for the `B` matrix is as follows:
     *      B_IAJB = (AI||BJ)
     *
     *  @param ghf_structure        The GHF QCStructure which contains the ground state parameters of the performed GHF calculation.
     *  @param gsq_hamiltonian      The generalised, second quantized hamiltonian, which contains the necessary two electron operators.
     */
    const GQCP::Matrix<Scalar> calculatePartialStabilityMatrixB(const QCStructure<QCModel::GHF<Scalar>>& ghf_structure, const GSQHamiltonian<Scalar>& gsq_hamiltonian) const {

        // Get the ground state parameters from the given QCStructure.
        const auto& parameters = ghf_structure.groundStateParameters();

        // Determine the number of occupied and virtual orbitals.
        const auto& number_of_occupied_orbitals = parameters.numberOfElectrons();
        const auto& number_of_virtual_orbitals = gsq_hamiltonian.numberOfOrbitals() - number_of_occupied_orbitals;

        // We need the two-electron integrals in MO basis, hence why we transform them with the coefficient matrix.
        // The ground state coefficient matrix is obtained from the QCModel.
        // We need the anti-symmetrized tensor: (AI||BJ) = (AI|BJ) - (AJ|BI). This is obtained by the `.antisymmetrized()` method.
        std::cout << parameters.coefficientMatrix() << std::endl;
        const auto& g = gsq_hamiltonian.twoElectron().transformed(parameters.coefficientMatrix()).antisymmetrized();

        // The next step is to create the needed tensor slice.
        GQCP::Tensor<Scalar, 4> B_iajb(number_of_occupied_orbitals, number_of_virtual_orbitals, number_of_occupied_orbitals, number_of_virtual_orbitals);
        for (int i = 0; i < number_of_occupied_orbitals; i++) {
            for (int a = number_of_occupied_orbitals; a < parameters.numberOfSpinors(); a++) {
                for (int j = 0; j < number_of_occupied_orbitals; j++) {
                    for (int b = number_of_occupied_orbitals; b < parameters.numberOfSpinors(); b++) {
                        B_iajb(i, a - number_of_occupied_orbitals, j, b - number_of_occupied_orbitals) = g.parameters()(a, i, b, j);  // -0.149 & 0.419 switched atm. tensor elements right place, rework reshape function.
                    }
                }
            }
        }
        B_iajb.print();
        //Finally, reshape the tensor to a matrix.
        const auto B_matrix = B_iajb.reshape(number_of_occupied_orbitals * number_of_virtual_orbitals, number_of_occupied_orbitals * number_of_virtual_orbitals);

        return B_matrix;
    }


    /**
     *  @return the internal stability matrix of the real GHF method.
     *
     *  @note The internal stability condition of the real GHF method is checked using A+B.
     *
     *  @param ghf_structure        The GHF QCStructure which contains the ground state parameters of the performed GHF calculation.
     *  @param gsq_hamiltonian      The generalised, second quantized hamiltonian, which contains the necessary two electron operators.
     */
    template <typename S = Scalar, typename = IsReal<S>>
    const GQCP::Matrix<double> internalStabilityMatrix(const QCStructure<QCModel::GHF<double>>& ghf_structure, const GSQHamiltonian<double>& gsq_hamiltonian) {
        return calculatePartialStabilityMatrixA(ghf_structure, gsq_hamiltonian) + calculatePartialStabilityMatrixB(ghf_structure, gsq_hamiltonian);
    }


    /**
     *  @return the external stability matrix of the real GHF method.
     *
     *  @note The external stability condition of the real GHF method is checked using A-B.
     *
     *  @param ghf_structure        The GHF QCStructure which contains the ground state parameters of the performed GHF calculation.
     *  @param gsq_hamiltonian      The generalised, second quantized hamiltonian, which contains the necessary two electron operators.
     */
    template <typename S = Scalar, typename = IsReal<S>>
    const GQCP::Matrix<double> externalStabilityMatrix(const QCStructure<QCModel::GHF<double>>& ghf_structure, const GSQHamiltonian<double>& gsq_hamiltonian) {
        return calculatePartialStabilityMatrixA(ghf_structure, gsq_hamiltonian) - calculatePartialStabilityMatrixB(ghf_structure, gsq_hamiltonian);
    }


    /**
     *  @return the internal stability matrix of the complex GHF method.
     *
     *  @note The internal stability condition of the real GHF method is checked using (A,   B  )
     *                                                                                 (B^*, A^*)
     *
     *  @param ghf_structure        The GHF QCStructure which contains the ground state parameters of the performed GHF calculation
     *  @param gsq_hamiltonian      The generalised, second quantized hamiltonian, which contains the necessary two electron operators.
     */
    template <typename S = Scalar, typename = IsComplex<S>>
    const GQCP::Matrix<complex> internalStabilityMatrix(const QCStructure<QCModel::GHF<complex>>& ghf_structure, const GSQHamiltonian<double>& gsq_hamiltonian) {

        // Calculate the necessary partial stability matrices.
        const auto A = calculatePartialStabilityMatrixA(ghf_structure, gsq_hamiltonian);
        const auto B = calculatePartialStabilityMatrixB(ghf_structure, gsq_hamiltonian);

        // Determine the dimensions of the total stability matrix.
        const auto K = A.dimension(0);
        const auto dim = 2 * K;

        // Create the total stability matrix as specified above in the documentation.
        auto H = GQCP::Matrix<Scalar> {dim, dim};

        H.topLeftCorner(K, K) = A;
        H.topRightCorner(K, K) = B;
        H.bottomLeftCorner(K, K) = B.conjugate();
        H.bottomRightCorner(K, K) = A.conjugate();

        return H;
    }
};

}  // namespace QCMethod
}  // namespace GQCP
