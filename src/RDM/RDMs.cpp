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
#include "RDM/RDMs.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  Constructor with @param one_rdm
 *  where one_rdm_aa and one_rdm_bb are @param one_rdm/2
 */
OneRDMs::OneRDMs(const OneRDM& one_rdm) :
        one_rdm (one_rdm),
        one_rdm_aa (one_rdm.get_matrix_representation()/2),
        one_rdm_bb (one_rdm.get_matrix_representation()/2)
{}

/**
 *  Constructor with @param one_rdm_aa and @param one_rdm_bb
 *  were one_rdm = @param one_rdm_aa + @param one_rdm_bb
 */
OneRDMs::OneRDMs(const OneRDM& one_rdm_aa, const OneRDM& one_rdm_bb) :
        one_rdm (one_rdm_aa.get_matrix_representation() + one_rdm_bb.get_matrix_representation()),
        one_rdm_aa (one_rdm_aa),
        one_rdm_bb (one_rdm_bb)
{}

/**
 *  Constructor
 *  where two_rdm = @param two_rdm_aaaa + @param two_rdm_aabb + @param two_rdm_bbaa  + @param two_rdm_bbbb
 */
TwoRDMs::TwoRDMs(const TwoRDM& two_rdm_aaaa, const TwoRDM& two_rdm_aabb, const TwoRDM& two_rdm_bbaa, const TwoRDM& two_rdm_bbbb) :
        two_rdm (two_rdm_aaaa.get_matrix_representation() + two_rdm_aabb.get_matrix_representation() + two_rdm_bbaa.get_matrix_representation() + two_rdm_bbbb.get_matrix_representation()),
        two_rdm_aaaa (two_rdm_aaaa),
        two_rdm_aabb (two_rdm_aabb),
        two_rdm_bbaa (two_rdm_bbaa),
        two_rdm_bbbb (two_rdm_bbbb)
{}


}  // namespace GQCP
