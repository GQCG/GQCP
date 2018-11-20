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
#ifndef GQCP_BASEHAMILTONIANPARAMETERS_HPP
#define GQCP_BASEHAMILTONIANPARAMETERS_HPP


#include <memory>

#include "AOBasis.hpp"


namespace GQCP {


/**
 *  A base class for representing Hamiltonian parameters, i.e. the one- and two-electron integrals in the second-quantized expression of the Hamiltonian
 */
class BaseHamiltonianParameters {
protected:
    double scalar;  // a scalar interaction term
    std::shared_ptr<GQCP::AOBasis> ao_basis;  // the initial atomic orbitals

public:
    // CONSTRUCTORS
    /**
     *  @param ao_basis     the initial AO basis
     *  @param scalar       the scalar interaction term
     */
    BaseHamiltonianParameters(std::shared_ptr<GQCP::AOBasis> ao_basis, double scalar=0.0);


    // DESTRUCTOR
    /**
     *  Provide a pure virtual destructor to make the class abstract
     */
    virtual ~BaseHamiltonianParameters() = 0;

    // GETTERS
    const std::shared_ptr<GQCP::AOBasis>& get_ao_basis() const { return this->ao_basis; }
};


}  // namespace GQCP


#endif  // GQCP_BASEHAMILTONIANPARAMETERS_HPP
