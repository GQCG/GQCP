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
#include "RDM/RDMs.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  A struct that holds the spin-summed, as well as the spin-resolved 1-RDMs
 */
OneRDMs::OneRDMs(const OneRDM& one_rdm) :
        one_rdm (one_rdm),
        one_rdm_aa (one_rdm.get_matrix_representation()/2),
        one_rdm_bb (one_rdm.get_matrix_representation()/2)
{}


/**
 *  A constructor that creates the total 1-RDM as the sum of the spin-resolved 1-RDMs
 *
 *  @param one_rdm_aa       the alpha-alpha 1-RDM
 *  @param one_rdm_bb       the beta-beta 1-RDM
 */
OneRDMs::OneRDMs(const OneRDM& one_rdm_aa, const OneRDM& one_rdm_bb) :
        one_rdm (one_rdm_aa.get_matrix_representation() + one_rdm_bb.get_matrix_representation()),
        one_rdm_aa (one_rdm_aa),
        one_rdm_bb (one_rdm_bb)
{}


/**
 *  A constructor that creates the total 2-RDM as the sum of the spin-resolved 2-RDMs
 *
 *  @param two_rdm_aaaa     the alpha-alpha-alpha-alpha 2-RDM
 *  @param two_rdm_aabb     the alpha-alpha-beta-beta 2-RDM
 *  @param two_rdm_bbaa     the beta-beta-alpha-alpha 2-RDM
 *  @param two_rdm_bbbb     the beta-beta-beta-beta 2-RDM
 */
TwoRDMs::TwoRDMs(const TwoRDM& two_rdm_aaaa, const TwoRDM& two_rdm_aabb, const TwoRDM& two_rdm_bbaa, const TwoRDM& two_rdm_bbbb) :
        two_rdm (two_rdm_aaaa.get_matrix_representation() + two_rdm_aabb.get_matrix_representation() + two_rdm_bbaa.get_matrix_representation() + two_rdm_bbbb.get_matrix_representation()),
        two_rdm_aaaa (two_rdm_aaaa),
        two_rdm_aabb (two_rdm_aabb),
        two_rdm_bbaa (two_rdm_bbaa),
        two_rdm_bbbb (two_rdm_bbbb)
{}


}  // namespace GQCP
