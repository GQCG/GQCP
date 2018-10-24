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
#ifndef GQCG_BASERDMBUILDER_HPP
#define GQCG_BASERDMBUILDER_HPP


#include "RDM/OneRDM.hpp"
#include "RDM/TwoRDM.hpp"
#include "RDM/RDMs.hpp"
#include "FockSpace/BaseFockSpace.hpp"


namespace GQCP {


/**
 *  BaseRDMBuilder is an abstract base class for the calculation of a density matrix from a given wave function
 *  or coefficient expansion in a Fock space
 */
class BaseRDMBuilder {
public:
    // CONSTRUCTOR
    BaseRDMBuilder() = default;


    // DESTRUCTOR
    /**
     *  Provide a pure virtual destructor to make the class abstract
     */
    virtual ~BaseRDMBuilder() = 0;


    // PURE VIRTUAL PUBLIC METHODS
    /**
     *  @return all 1-RDMs from a coefficient vector @param x
     */
    virtual OneRDMs calculate1RDMs(const Eigen::VectorXd& x) = 0;

    /**
     *  @return all 2-RDMs from a coefficient vector @param x
     */
    virtual TwoRDMs calculate2RDMs(const Eigen::VectorXd& x) = 0;

    /**
     *  @return the Fock space of the RDMBuilder
     */
    virtual BaseFockSpace* get_fock_space() = 0;
};


}  // namespace GQCG


#endif  // GQCG_BASERDMBUILDER_HPP
