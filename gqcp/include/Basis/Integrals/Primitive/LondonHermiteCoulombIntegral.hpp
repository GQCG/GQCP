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


#include "Mathematical/Representation/Matrix.hpp"
#include "Utilities/complex.hpp"


namespace GQCP {


/**
 *  An implementation of the (auxiliary) London Hermite Coulomb integral R^{k1, n}_{tuv}.
 */
class LondonHermiteCoulombIntegral {
private:
    // The exponent of the Hermite Gaussian.
    double p;

    // The k-vector of the London overlap distribution.
    Vector<double, 3> k1;

    // The center of the Hermite Gaussian.
    Vector<double, 3> P;

    // The center of the Coulomb potential.
    Vector<double, 3> C;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  @param k1           The k-vector of the London overlap distribution.
     *  @param p            The exponent of the Hermite Gaussian.
     *  @param P            The center of the Hermite Gaussian.
     *  @param C            The center of the Coulomb potential.
     */
    LondonHermiteCoulombIntegral(const Vector<double, 3>& k1, const double p, const Vector<double, 3>& P, const Vector<double, 3>& C);


    /*
     *  MARK: London Hermite Coulomb integral implementation
     */

    /**
     *  Calculate the value for the (auxiliary) London Hermite Coulomb integral R^{k1, n}_{tuv}(p, P, C).
     * 
     *  @param n            The order of the London Hermite Coulomb integral, i.e. the order of the Boys function.
     *  @param t            The derivative degree in P_x.
     *  @param u            The derivative degree in P_y.
     *  @param v            The derivative degree in P_z.
     */
    complex operator()(const size_t n, const int t, const int u, const int v) const;
};


}  // namespace GQCP
