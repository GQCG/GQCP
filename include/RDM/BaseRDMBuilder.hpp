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
 *  BaseRDMBuilder is an abstract base class whose derived classes are capable of calculating 1- and 2-RDMs from wave functions
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


    // PURE VIRTUAL GETTERS
    virtual BaseFockSpace* get_fock_space() = 0;


    // PURE VIRTUAL PUBLIC METHODS
    /**
     *  @param x        the coefficient vector representing the wave function
     *
     *  @return all 1-RDMs given a coefficient vector
     */
    virtual OneRDMs calculate1RDMs(const Eigen::VectorXd& x) = 0;

    /**
     *  @param x        the coefficient vector representing the wave function
     *
     *  @return all 2-RDMs given a coefficient vector
     */
    virtual TwoRDMs calculate2RDMs(const Eigen::VectorXd& x) = 0;
};


}  // namespace GQCG


#endif  // GQCG_BASERDMBUILDER_HPP
