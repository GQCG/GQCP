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




#include "Basis/SpinorBasis/JacobiRotationParameters.hpp"
#include "Basis/TransformationMatrix.hpp"
#include "Mathematical/Representation/QCMatrix.hpp"
#include "Mathematical/ScalarFunction.hpp"
#include "Processing/RDM/OneRDM.hpp"
#include "Processing/RDM/TwoRDM.hpp"
#include "Utilities/type_traits.hpp"
#include "Basis/SpinorBasis/SpinComponent.hpp"

#include <array>


namespace GQCP {


/**
 *  A class that represents an unrestricted second-quantized one-electron operator: it holds the matrix representation of its parameters for both spin components, which are (usually) integrals over first-quantized operators
 *
 *  @tparam _Scalar             the scalar type, i.e. the scalar representation of one of the parameters
 *  @tparam _Components         the number of components of the second-quantized operator
 */
template <typename _Scalar, size_t _Components>
class USQOneElectronOperator{
public:

    using Scalar = _Scalar;
    static constexpr auto Components = _Components;

private:
    std::array<QCMatrix<Scalar>, Components> fs;  // all the matrix representations (hence the s) of the parameters (integrals) of the different components of this second-quantized operator

public:

    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param fs_alpha    all the matrix representations of the alpha spin component (hence the s) of the parameters (integrals) of the different components of this second-quantized operator
     *  @param fs_beta     all the matrix representations of the beta spin component (hence the s) of the parameters (integrals) of the different components of this second-quantized operator
     */
    USQOneElectronOperator(const std::array<QCMatrix<Scalar>, Components>& fs_alpha, const std::array<QCMatrix<Scalar>, Components>& fs_beta) : 
        fs_alpha (fs_alpha),
        fs_beta (fs_beta)
    {
        // Check if the given matrix representations have the same dimensions
        const auto dimension_of_first_alpha = this->fs_alpha[0].dimension();
        const auto dimension_of_first_beta = this->fs_beta[0].dimension();
        for (size_t i = 1; i < Components; i++) {

            const auto dimension_of_ith_alpha = this->fs_alpha[i].dimension();
            const auto dimension_of_ith_beta = this->fs_beta[i].dimension();

            if (dimension_of_first_alpha != dimension_of_ith_alpha || dimension_of_first_beta != dimension_of_ith_beta){
                throw std::invalid_argument("USQOneElectronOperator(const std::array<QCMatrix<Scalar>, Components>&, const std::array<QCMatrix<Scalar>, components>&): 
                The given matrix representations do not have the same dimensions for either the alpha or beta component.");
            }
        }
    }    


    /**
     *  A constructor for ScalarUSQOneElectronOperators that doesn't require the argument to be a vector of just one element.
     * 
     *  @param f_alpha            the matrix representation of the integrals of this scalar second-quantized operator, for the alpha spin component
     *  @param f_beta             the matrix representation of the integrals of this scalar second-quantized operator, for the beta spin component
     * 
     *  @note This constructor is only available for ScalarUSQOneElectronOperators (for the std::enable_if, see https://stackoverflow.com/a/17842695/7930415)
     */
    template <size_t Z = Components>
    USQOneElectronOperator(const QCMatrix<Scalar>& f_alpha, typename std::enable_if<Z == 1>::type* = 0, const QCMatrix<Scalar>& f_beta, typename std::enable_if<Z == 1>::type* = 0) :
        USQOneElectronOperator(std::array<QCMatrix<Scalar>, 1>{f_alpha}, std::array<QCMatrix<Scalar>, 1>{f_beta})
    {}


    /**
     *  Construct a one-electron operator with parameters that are zero, for both the alpha and beta spin components
     * 
     *  @param dim          the dimension of the matrix representation of the parameters, i.e. the number of orbitals/sites. This dimension is the same for the alpha and beta component.
     */
    USQOneElectronOperator(const size_t dim, const size_t dim) {
        for (size_t i = 0; i < Components; i++) {
            this->fs_alpha[i] = QCMatrix<Scalar>::Zero(dim, dim);
            this->fs_beta[i] = QCMatrix<Scalar>::Zero(dim, dim);
        }
    }


    /**
     *  Default constructor: construct a one-electron operator with parameters that are zero, for both spin components.
     */
    USQOneElectronOperator() :
        USQOneElectronOperator(0, 0)  // dimensions of the representations are zero
    {}


    /*
     *  GETTERS
     */
    size_t get_dim_alpha() const { return this->f_alpha.dimension(); }
    size_t get_dim_beta() const { return this->f_beta.dimension(); }
    size_t get_K_alpha() const { return this->f_alpha.dimension(); }
    size_t get_K_beta() const { return this->f_beta.dimension(); }


    /*
     *  OPERATORS
     */
    
}

} // namespace GQCP