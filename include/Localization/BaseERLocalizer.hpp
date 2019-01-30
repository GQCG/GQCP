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
#ifndef BaseERLocalizer_hpp
#define BaseERLocalizer_hpp


#include "HamiltonianParameters/HamiltonianParameters.hpp"


namespace GQCP {


/**
 *  A base class to implement the Edmiston-Ruedenberg localization method
 */
class BaseERLocalizer {
protected:
    size_t N_P;  // the number of electron pairs
    double threshold;  // the threshold for maximization on subsequent localization indices
    size_t maximum_number_of_iterations;  // the maximum number of iterations for the localization algorithm

    bool is_converged = false;
    size_t iterations = 0;  // the number of iterations

public:
    // CONSTRUCTORS
    /**
     *  @param N_P                              the number of electron pairs
     *  @param threshold                        the threshold for maximization on subsequent localization indices
     *  @param maximum_number_of_iterations     the maximum number of iterations for the localization algorithm
     */
    BaseERLocalizer(size_t N_P, double threshold=1.0e-08, size_t maximum_number_of_iterations=128);


    // PUBLIC METHODS
    /**
     *  Localize the Hamiltonian parameters by maximizing the Edmiston-Ruedenberg localization index
     *
     *  @param ham_par      the Hamiltonian parameters (in an orthonormal basis) that should be localized
     */
    virtual void localize(HamiltonianParameters& ham_par) = 0;
};



}  // namespace GQCP


#endif /* BaseERLocalizer_hpp */
