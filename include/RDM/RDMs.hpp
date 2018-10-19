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


struct OneRDMs {

    OneRDM one_rdm;  // spin-summed (total) 1-RDM

    OneRDM one_rdm_aa;  // alpha-alpha (a-a) 1-RDM
    OneRDM one_rdm_bb;  // beta-beta (b-b) 1-RDM

    // CONSTRUCTORS
    /**
     *  Constructor with @param one_rdm
     *  where one_rdm_aa and one_rdm_bb are @param one_rdm/2
     */
    OneRDMs(const OneRDM& one_rdm);

    /**
     *  Constructor with @param one_rdm_aa and @param one_rdm_bb
     *  were one_rdm = @param one_rdm_aa + @param one_rdm_bb
     */
    OneRDMs(const OneRDM& one_rdm_aa, const OneRDM& one_rdm_bb);

};


struct TwoRDMs {

    TwoRDM two_rdm;  // spin-summed (total) 2-RDM

    TwoRDM two_rdm_aaaa;  // a-a-a-a 2-RDM
    TwoRDM two_rdm_aabb;  // a-a-b-b 2-RDM
    TwoRDM two_rdm_bbaa;  // b-a-a-b 2-RDM
    TwoRDM two_rdm_bbbb;  // b-b-b-b 2-RDM


    // CONSTRUCTORS
    /**
     *  Constructor
     *  where two_rdm = @param two_rdm_aaaa + @param two_rdm_aabb + @param two_rdm_bbaa  + @param two_rdm_bbbb
     */
    TwoRDMs(const TwoRDM& two_rdm_aaaa, const TwoRDM& two_rdm_aabb, const TwoRDM& two_rdm_bbaa, const TwoRDM& two_rdm_bbbb);

};


}  // namespace GQCP

#endif  // GQCP_RDMS_HPP
