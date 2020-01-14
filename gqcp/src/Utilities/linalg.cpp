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
#include "Utilities/linalg.hpp"


namespace GQCP {


/**
 *  @param eigenvalues1     the first set of eigenvalues
 *  @param eigenvalues2     the second set of eigenvalues
 *  @param tolerance        the tolerance for comparison
 *
 *  @return if two sets of eigenvalues are equal within a given tolerance
 */
bool areEqualEigenvalues(const VectorX<double>& eigenvalues1, const VectorX<double>& eigenvalues2, double tolerance) {
    return eigenvalues1.isApprox(eigenvalues2, tolerance);
}


/**
 *  @param eigenvector1     the first eigenvector
 *  @param eigenvector2     the second eigenvector
 *  @param tolerance        the tolerance for comparison
 *
 *  @return if two eigenvectors are equal within a given tolerance
 */
bool areEqualEigenvectors(const VectorX<double>& eigenvector1, const VectorX<double>& eigenvector2, double tolerance) {

    //  Eigenvectors are equal if they are equal up to their sign.
    return (eigenvector1.isApprox(eigenvector2, tolerance) || eigenvector1.isApprox(-eigenvector2, tolerance));
}


/**
 *  @param eigenvectors1        the first set of eigenvectors
 *  @param eigenvectors2        the second set of eigenvectors
 *  @param tolerance            the tolerance for comparison
 *
 *  @return if two sets of eigenvectors are equal within a given tolerance
 */
bool areEqualSetsOfEigenvectors(const MatrixX<double>& eigenvectors1, const MatrixX<double>& eigenvectors2, double tolerance) {

    // Check if the dimensions of the eigenvectors are equal.
    if (eigenvectors1.cols() != eigenvectors2.cols()) {
        throw std::invalid_argument("areEqualSetsOfEigenvectors(MatrixX<double>, MatrixX<double>, double): Cannot compare the two sets of eigenvectors as they have different dimensions.");
    }

    if (eigenvectors1.rows() != eigenvectors2.rows()) {
        throw std::invalid_argument("areEqualSetsOfEigenvectors(MatrixX<double>, MatrixX<double>, double): Cannot compare the two sets of eigenvectors as they have different dimensions.");
    }


    for (size_t i = 0; i < eigenvectors1.cols(); i++) {
        const VectorX<double> eigenvector1 = eigenvectors1.col(i);
        const VectorX<double> eigenvector2 = eigenvectors2.col(i);

        if (!areEqualEigenvectors(eigenvector1, eigenvector2, tolerance)) {
            return false;
        }
    }

    return true;
}


}  // namespace GQCP
