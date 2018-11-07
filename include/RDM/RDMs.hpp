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
#ifndef GQCP_RDMS_HPP
#define GQCP_RDMS_HPP


#include "RDM/OneRDM.hpp"
#include "RDM/TwoRDM.hpp"


namespace GQCP {


/**
 *  A struct that holds the spin-summed, as well as the spin-resolved 1-RDMs
 */
struct OneRDMs {
    OneRDM one_rdm;  // spin-summed (total) 1-RDM

    OneRDM one_rdm_aa;  // alpha-alpha (a-a) 1-RDM
    OneRDM one_rdm_bb;  // beta-beta (b-b) 1-RDM


    // CONSTRUCTORS
    /**
     *  A constructor that creates the spin-resolved 1-RDMs as half of the total 1-RDM
     *
     *  @param one_rdm      the spin-summed 1-RDM
     */
    OneRDMs(const OneRDM& one_rdm);

    /**
     *  A constructor that creates the total 1-RDM as the sum of the spin-resolved 1-RDMs
     *
     *  @param one_rdm_aa       the alpha-alpha 1-RDM
     *  @param one_rdm_bb       the beta-beta 1-RDM
     */
    OneRDMs(const OneRDM& one_rdm_aa, const OneRDM& one_rdm_bb);

};


/**
 *  A struct that holds the spin-summed, as well as the spin-separated 2-RDMs
 */
struct TwoRDMs {
    TwoRDM two_rdm;  // spin-summed (total) 2-RDM

    TwoRDM two_rdm_aaaa;  // a-a-a-a 2-RDM
    TwoRDM two_rdm_aabb;  // a-a-b-b 2-RDM
    TwoRDM two_rdm_bbaa;  // b-a-a-b 2-RDM
    TwoRDM two_rdm_bbbb;  // b-b-b-b 2-RDM


    // CONSTRUCTORS
    /**
     *  A constructor that creates the total 2-RDM as the sum of the spin-resolved 2-RDMs
     *
     *  @param two_rdm_aaaa     the alpha-alpha-alpha-alpha 2-RDM
     *  @param two_rdm_aabb     the alpha-alpha-beta-beta 2-RDM
     *  @param two_rdm_bbaa     the beta-beta-alpha-alpha 2-RDM
     *  @param two_rdm_bbbb     the beta-beta-beta-beta 2-RDM
     */
    TwoRDMs(const TwoRDM& two_rdm_aaaa, const TwoRDM& two_rdm_aabb, const TwoRDM& two_rdm_bbaa, const TwoRDM& two_rdm_bbbb);

};


}  // namespace GQCP

#endif  // GQCP_RDMS_HPP
