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
#include "WaveFunction/WaveFunction.hpp"
#include "FockSpace/ProductFockSpace.hpp"

namespace GQCP {


/*
 * CONSTRUCTORS
 */

/**
 *  Construct a normalized wave function from possibly non-normalized coefficients
 *
 *  @param base_fock_space      the Fock space in which the wave function 'lives'
 *  @param coefficients         the expansion coefficients
 */
WaveFunction::WaveFunction(const BaseFockSpace& base_fock_space, const VectorX<double>& coefficients) :
    fock_space (BaseFockSpace::CloneToHeap(base_fock_space)),
    coefficients (coefficients)
{
    if (std::abs(this->coefficients.norm() - 1.0) > 1.0e-12) {  // normalize the coefficients if they aren't
        this->coefficients.normalize();
    }
}




/*
 *  PUBLIC METHODS
 */
/**
 *  @return the Shannon entropy (or information content) of the wave function
 */
double WaveFunction::calculateShannonEntropy() const {

    // Sum over the Fock space dimension, and only include the term if c_k != 0
    // We might as well replace all coeffients that are 0 by 1, since log(1) = 0 so there is no influence on the final entropy value
    Eigen::ArrayXd coefficients_replaced = this->coefficients.unaryExpr([](double c) { return c < 1.0e-18 ? 1 : c;});  // replace 0 by 1

    Eigen::ArrayXd coefficients_squared = coefficients_replaced.square();
    Eigen::ArrayXd log_coefficients_squared = coefficients_squared.log();  // natural logarithm (ln)

    return - 1 / std::log(2) * (coefficients_squared * log_coefficients_squared).sum();
}


/**
 *  Transform the underlying ONV basis of the wave function (only for FCI [ProductFockSpace])
 *
 *  @param T    the transformation matrix between the old and the new orbital basis
 */
void WaveFunction::basisTransform(const SquareMatrix<double>& T) {

    if (fock_space->get_type() != FockSpaceType::ProductFockSpace) {
        throw std::invalid_argument("WaveFunction::basisTransform(SquareMatrix<double>): This is not an FCI wave function");
    }

    auto K = fock_space->get_K();

    if (K != T.get_dim()) {
        throw std::invalid_argument("WaveFunction::basisTransform(SquareMatrix<double>): number of orbitals does not match the dimension of the transformation matrix T");
    }

    Eigen::FullPivLU<Eigen::MatrixXd> LU_decomposer (T);

    SquareMatrix<double> L = SquareMatrix<double>::Zero(K, K);
    L.triangularView<Eigen::StrictlyLower>() = LU_decomposer.matrixLU();

    std::cout<<std::endl<<LU_decomposer.matrixLU()<<std::endl;
    std::cout<<std::endl<<L<<std::endl;

    SquareMatrix<double> U = SquareMatrix<double>(LU_decomposer.matrixLU().triangularView<Eigen::Upper>());

    std::cout<<std::endl<<U<<std::endl;
    std::cout<<std::endl<<U.inverse()<<std::endl;

    SquareMatrix<double> t =  - L + U.inverse();

    const auto& product_fock_space = dynamic_cast<const ProductFockSpace&>(*fock_space);
    const FockSpace& fock_space_alpha = product_fock_space.get_fock_space_alpha();
    const FockSpace& fock_space_beta = product_fock_space.get_fock_space_beta();

    auto dim_alpha = fock_space_alpha.get_dimension();
    auto dim_beta = fock_space_beta.get_dimension();
    auto N_alpha = fock_space_alpha.get_N();
    auto N_beta = fock_space_beta.get_N();

    VectorX<double> current_coefficients = this->coefficients;
    VectorX<double> correction_coefficients = VectorX<double>::Zero(product_fock_space.get_dimension());

    std::cout<<std::endl<<"t:"<<t<<std::endl;

    for (size_t m = 0; m < K; m++) {

        ONV alpha = fock_space_alpha.makeONV(0);

        for (size_t I_alpha = 0; I_alpha < dim_alpha; I_alpha++) {
            if (!alpha.isOccupied(m)) {
                for (size_t e1 = 0; e1 < N_alpha; e1++) {  // e1 (electron 1) loops over the (number of) electrons
                    size_t p = alpha.get_occupation_index(e1);  // retrieve the index of a given electron

                    if (p < m) {

                        size_t address = I_alpha - fock_space_alpha.get_vertex_weights(p, e1 + 1);
                        size_t e2 = e1 + 1;
                        size_t q = p + 1;
                        int sign = 1;

                        fock_space_alpha.shiftUntilNextUnoccupiedOrbital<1>(alpha, address, q, e2, sign);
                        while (q != m) {
                            q++;
                            fock_space_alpha.shiftUntilNextUnoccupiedOrbital<1>(alpha, address, q, e2, sign);
                        }

                        address += fock_space_alpha.get_vertex_weights(q, e2);

                        for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {
                            correction_coefficients(I_alpha * dim_beta + I_beta) +=
                                    sign * t(p, m) * current_coefficients(address * dim_beta + I_beta);
                        }
                    }

                    if (p > m) {

                        size_t address = I_alpha - fock_space_alpha.get_vertex_weights(p, e1 + 1);
                        size_t e2 = e1 - 1;
                        size_t q = p - 1;
                        int sign = 1;

                        fock_space_alpha.shiftUntilPreviousUnoccupiedOrbital<1>(alpha, address, q, e2, sign);
                        while (q != m) {
                            q--;
                            fock_space_alpha.shiftUntilPreviousUnoccupiedOrbital<1>(alpha, address, q, e2, sign);
                        }

                        address += fock_space_alpha.get_vertex_weights(q, e2 + 2);

                        std::cout<<I_alpha<<" : "<<"p :"<<p<<" m :"<<m<<" address: "<<address<<" sign: "<<sign<<std::endl;
                        for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {
                            correction_coefficients(I_alpha * dim_beta + I_beta) +=
                                    sign * t(p, m) * current_coefficients(address * dim_beta + I_beta);
                        }
                    }
                }


            } else {
                for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {
                    correction_coefficients(I_alpha * dim_beta + I_beta) +=
                            (t(m, m) - 1) * current_coefficients(I_alpha * dim_beta + I_beta);
                }
            }


            if (I_alpha < dim_alpha - 1) {  // prevent the last permutation to occur
                fock_space_alpha.setNextONV(alpha);
            }
        }


        current_coefficients += correction_coefficients;

        std::cout<<std::endl<<"alpha "<<m<<" :"<<current_coefficients<<std::endl;
        correction_coefficients.setZero();

        // BETA-BRANCH

        ONV beta = fock_space_beta.makeONV(0);

        for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {
            if (!beta.isOccupied(m)) {
                for (size_t e1 = 0; e1 < N_beta; e1++) {  // e1 (electron 1) loops over the (number of) electrons
                    size_t p = beta.get_occupation_index(e1);  // retrieve the index of a given electron

                    if (p < m) {

                        size_t address = I_beta - fock_space_beta.get_vertex_weights(p, e1 + 1);
                        size_t e2 = e1 + 1;
                        size_t q = p + 1;
                        int sign = 1;

                        fock_space_beta.shiftUntilNextUnoccupiedOrbital<1>(beta, address, q, e2, sign);
                        while (q != m) {
                            q++;
                            fock_space_beta.shiftUntilNextUnoccupiedOrbital<1>(beta, address, q, e2, sign);
                        }

                        address += fock_space_beta.get_vertex_weights(q, e2);

                        for (size_t I_alpha = 0; I_alpha < dim_alpha; I_alpha++) {
                            correction_coefficients(I_alpha * dim_beta + I_beta) +=
                                    sign * t(p, m) * current_coefficients(I_alpha * dim_beta + address);
                        }
                    }

                    if (p > m) {

                        size_t address = I_beta - fock_space_beta.get_vertex_weights(p, e1 + 1);
                        size_t e2 = e1 - 1;
                        size_t q = p - 1;

                        int sign = 1;

                        fock_space_beta.shiftUntilPreviousUnoccupiedOrbital<1>(beta, address, q, e2, sign);
                        while (q != m) {
                            q--;
                            fock_space_beta.shiftUntilPreviousUnoccupiedOrbital<1>(beta, address, q, e2, sign);
                        }

                        address += fock_space_beta.get_vertex_weights(q, e2 + 2);

                        for (size_t I_alpha = 0; I_alpha < dim_alpha; I_alpha++) {
                            correction_coefficients(I_alpha * dim_beta + I_beta) +=
                                    sign * t(p, m) * current_coefficients(I_alpha * dim_beta + address);
                        }
                    }
                }
            } else {
                for (size_t I_alpha = 0; I_alpha < dim_alpha; I_alpha++) {
                    correction_coefficients(I_alpha * dim_beta + I_beta) +=
                            (t(m, m) - 1) * current_coefficients(I_alpha * dim_beta + I_beta);
                }
            }

            if (I_beta < dim_beta - 1) {  // prevent the last permutation to occur
                fock_space_beta.setNextONV(beta);
            }
        }

        current_coefficients += correction_coefficients;
        correction_coefficients.setZero();
    }

    this->coefficients = current_coefficients;
}

}  // namespace GQCP
