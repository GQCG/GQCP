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


#include "Mathematical/Representation/SquareMatrix.hpp"
#include "Operator/SecondQuantized/ModelHamiltonian/AdjacencyMatrix.hpp"


namespace GQCP {


/**
 *  The Hubbard hopping matrix.
 *
 *  @tparam _Scalar         The scalar type of the elements of the hopping matrix.
 */
template <typename _Scalar>
class HoppingMatrix {
public:
    // The scalar type of the elements of the hopping matrix.
    using Scalar = _Scalar;


private:
    // The matrix representation of the hopping matrix.
    SquareMatrix<Scalar> H;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  Create a hopping matrix from its representation as a `SquareMatrix`.
     *
     *  @param H        The Hubbard hopping matrix, represented as a `SquareMatrix`.
     */
    HoppingMatrix(const SquareMatrix<Scalar>& H) :
        H {H} {

        if (!H.isHermitian()) {
            throw std::invalid_argument("HoppingMatrix::HoppingMatrix(const SquareMatrix<Scalar>&): The given hopping matrix must be Hermitian.");
        }
    }


    /*
     *  MARK: Named constructors
     */

    /**
     *  Create a hopping matrix from a vector containing the upper triangle.
     *
     *  @param csline           A vector that contains the upper triangle (in column-major ordering) of the Hubbard hopping matrix.
     *
     *  @return The hopping matrix that corresponds to the given comma-separated line.
     */
    template <typename Z = Scalar>
    static enable_if_t<std::is_same<Z, double>::value, HoppingMatrix<double>> Dense(std::vector<double>& triagonal_data) {

        // Map the std::vector<double> into a VectorX<double> to be used into an other constructor
        GQCP::VectorX<double> upper_triangle = Eigen::Map<Eigen::VectorXd>(triagonal_data.data(), triagonal_data.size());
        return HoppingMatrix {SquareMatrix<double>::SymmetricFromUpperTriangleWithoutDiagonal(upper_triangle)};
    }


    /**
     *  Generate the Hubbard hopping matrix from an adjacency matrix and Hubbard model parameter t.
     *
     *  @param A        The adjacency matrix specifying the connectivity of the Hubbard lattice.
     *  @param t        The Hubbard parameter t. Note that a positive value for t means a negative neighbour hopping term.
     *
     *  @note This constructor is only available in the real case (for the std::enable_if, see https://stackoverflow.com/a/17842695/7930415).
     */
    template <typename Z = Scalar>
    static enable_if_t<std::is_same<Z, double>::value, HoppingMatrix<double>> Homogeneous(const AdjacencyMatrix& A, const double t) {
        return HoppingMatrix {-1 * t * A.matrix().cast<double>()};
    }


    /**
     *  Generate the Hubbard hopping matrix from an adjacency matrix and link vector containing values for different parameters t.
     *
     *  @param A                     The adjacency matrix specifying the connectivity of the Hubbard lattice.
     *  @param link_vector           A vector that contains the values of the links, in a given adjacency matrix.
     *
     *  @note This constructor is only available in the real case (for the std::enable_if, see https://stackoverflow.com/a/17842695/7930415).
     */
    template <typename Z = Scalar>
    static enable_if_t<std::is_same<Z, double>::value, HoppingMatrix<double>> FromLinkVector(const AdjacencyMatrix& A, std::vector<double>& link_vector) {
        // Generate the hopping (raw)matrix.
        size_t dim = A.matrix().rows();
        SquareMatrix<Scalar> H = SquareMatrix<Scalar>::Zero(dim);

        // Check if the dimension of the link vector matches the adjacency matrix.
        // If the adjacency matrix corressponds to a linear system.
        if (A.matrix()(0, dim - 1) == 0) {
            if (link_vector.size() != dim - 1) {
                throw std::invalid_argument("HoppingMatrix::FromLinkVector(const AdjacencyMatrix& A, std::vector<double>& link_vector): The dimension of the link vector does not match the adjacency matrix.");
            } else {
                for (size_t i = 0; i < dim - 1; i++) {  // row index
                    H(i, i + 1) = link_vector[i];
                    H(i + 1, i) = link_vector[i];
                }
            }
        }
        // If the adjacency matrix corresponds to a cyclic system.
        else if (A.matrix()(0, dim - 1) == 1) {
            if (link_vector.size() != dim) {
                throw std::invalid_argument("HoppingMatrix::FromLinkVector(const AdjacencyMatrix& A, std::vector<double>& link_vector): The dimension of the link vector does not match the adjacency matrix.");
            } else {
                for (size_t i = 0; i < dim - 1; i++) {  // row index
                    H(i, i + 1) = link_vector[i];
                    H(i + 1, i) = link_vector[i];
                }

                H(0, dim - 1) = link_vector[dim - 1];
                H(dim - 1, 0) = link_vector[dim - 1];
            }
        }

        return HoppingMatrix {H};
    }


    /*
     *  MARK: General information
     */

    /**
     *  @return The number of lattice sites corresponding used in this hopping matrix.
     */
    size_t numberOfLatticeSites() const { return this->matrix().dimension(); }


    /*
     *  MARK: Access
     */

    /**
     *  @return A read-only reference to the matrix representation of this hopping matrix.
     */
    const SquareMatrix<Scalar>& matrix() const { return this->H; }

    /**
     *  @return A writable reference to the matrix representation of this hopping matrix.
     */
    SquareMatrix<Scalar>& matrix() { return this->H; }
};


}  // namespace GQCP
