// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2018  the GQCG developers
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
#ifndef NRDMCalculator_hpp
#define NRDMCalculator_hpp


#include "FockSpace/FockSpace.hpp"

#include <iostream>
#include <algorithm>

#include <boost/range/adaptor/strided.hpp>
#include <boost/range/adaptor/sliced.hpp>
#include <boost/range/algorithm_ext/push_back.hpp>


namespace GQCP {


/**
 *  A class that can be used to calculate elements of N-th order density matrices for full CI wave functions (i.e. a FockSpace with spin orbitals)
 */
class NRDMCalculator {
private:
    FockSpace fock_space;  // the FCI Fock space with spin orbitals
    Eigen::VectorXd coeff;  // the expansion coefficient vector


public:
    // CONSTRUCTORS
    /**
     *  @param fock_space       the FCI Fock space with spin orbitals
     *  @param coeff            the expansion coefficient vector
     */
    NRDMCalculator(const FockSpace& fock_space, const Eigen::VectorXd& coeff);


    // OPERATORS
    /**
     *  @param indices_pack      the indices that specify the element of the N-RDM that has to be calculated
     */
    template<typename... size_ts>
    double operator()(size_ts... indices_pack) const {

        // Assume the user has given size_ts
        std::vector<size_t> indices {static_cast<size_t>(indices_pack)...};  // convert the pack to a vector so we can easily traverse


        if ((indices.size() == 0)) {
            return 1.0;  // assume the wave function is normalized
        }

        if ((indices.size() % 2) != 0) {
            throw std::invalid_argument("There must be an even number of indices as arguments.");
        }


        // Split the vector in ket (even) and bra (odd) indices
        std::vector<size_t> bra_indices;  // even
        boost::push_back(bra_indices, indices | boost::adaptors::strided(2));

        std::vector<size_t> ket_indices;  // odd
        boost::push_back(ket_indices, indices | boost::adaptors::sliced(1, indices.size()) | boost::adaptors::strided(2));
        std::reverse(ket_indices.begin(), ket_indices.end());


        return this->calculateElement(bra_indices, ket_indices);
    }


    // PUBLIC METHODS
    /**
     *  @param bra_indices      the indices of the orbitals that should be annihilated on the left (on the bra)
     *  @param ket_indices      the indices of the orbitals that should be annihilated on the right (on the ket)
     *
     *  @return an element of the N-RDM, as specified by the given bra and ket indices
     *
     *      calculateElement({0, 1}, {2, 1}) would calculate d^{(2)} (0, 1, 1, 2): the operator string would be a^\dagger_0 a^\dagger_1 a_2 a_1
     */
    double calculateElement(const std::vector<size_t>& bra_indices, const std::vector<size_t>& ket_indices) const;
};



}  // namespace GQCP

#endif /* NRDMCalculator_hpp */
