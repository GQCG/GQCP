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


namespace GQCP {


/**
 *  The spin-resolved two-electron density matrices.
 *
 *  @tparam Scalar      the scalar type
 */
template <typename Scalar>
struct SpinResolvedTwoDM {
    TwoDM<Scalar> two_rdm;  // spin-summed (total) 2-RDM

    TwoDM<Scalar> two_rdm_aaaa;  // a-a-a-a 2-RDM
    TwoDM<Scalar> two_rdm_aabb;  // a-a-b-b 2-RDM
    TwoDM<Scalar> two_rdm_bbaa;  // b-a-a-b 2-RDM
    TwoDM<Scalar> two_rdm_bbbb;  // b-b-b-b 2-RDM


    /*
     *  CONSTRUCTORS
     */

    /**
     *  A constructor that creates the total 2-RDM as the sum of the spin-resolved 2-RDMs
     *
     *  @param two_rdm_aaaa     the alpha-alpha-alpha-alpha 2-RDM
     *  @param two_rdm_aabb     the alpha-alpha-beta-beta 2-RDM
     *  @param two_rdm_bbaa     the beta-beta-alpha-alpha 2-RDM
     *  @param two_rdm_bbbb     the beta-beta-beta-beta 2-RDM
     */
    SpinResolvedTwoDM(const TwoDM<Scalar>& two_rdm_aaaa, const TwoDM<Scalar>& two_rdm_aabb, const TwoDM<Scalar>& two_rdm_bbaa, const TwoDM<Scalar>& two_rdm_bbbb) :
        two_rdm {two_rdm_aaaa.Eigen() + two_rdm_aabb.Eigen() + two_rdm_bbaa.Eigen() + two_rdm_bbbb.Eigen()},
        two_rdm_aaaa {two_rdm_aaaa},
        two_rdm_aabb {two_rdm_aabb},
        two_rdm_bbaa {two_rdm_bbaa},
        two_rdm_bbbb {two_rdm_bbbb} {}


    /**
     *  @return the dimension of the matrix representation of the 2-RDMs, i.e. the number of orbitals/sites
     */
    size_t dimension() const { return two_rdm.dimension(); }
};


}  // namespace GQCP
