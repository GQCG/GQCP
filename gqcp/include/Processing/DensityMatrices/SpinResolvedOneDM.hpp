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


#include "Processing/DensityMatrices/OneDM.hpp"


namespace GQCP {


/**
 *  The spin-resolved one-electron density matrices.
 *
 *  @tparam Scalar      the scalar type
 */
template <typename Scalar>
class SpinResolvedOneDM {
    OneDM<Scalar> one_rdm;  // spin-summed (total) 1-RDM

    OneDM<Scalar> one_rdm_aa;  // alpha-alpha (a-a) 1-RDM
    OneDM<Scalar> one_rdm_bb;  // beta-beta (b-b) 1-RDM


    /*
     *  CONSTRUCTORS
     */

    /**
     *  A constructor that creates the spin-resolved 1-RDMs as half of the total 1-RDM
     *
     *  @param one_rdm      the spin-summed 1-RDM
     */
    SpinResolvedOneDM(const OneDM<Scalar>& one_rdm) :
        one_rdm {one_rdm},
        one_rdm_aa {one_rdm / 2},
        one_rdm_bb {one_rdm / 2} {}


    /**
     *  A constructor that creates the total 1-RDM as the sum of the spin-resolved 1-RDMs
     *
     *  @param one_rdm_aa       the alpha-alpha 1-RDM
     *  @param one_rdm_bb       the beta-beta 1-RDM
     */
    SpinResolvedOneDM(const OneDM<Scalar>& one_rdm_aa, const OneDM<Scalar>& one_rdm_bb) :
        one_rdm {OneDM<double>(one_rdm_aa + one_rdm_bb)},
        one_rdm_aa {one_rdm_aa},
        one_rdm_bb {one_rdm_bb} {}


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return the dimension of the matrix representation of the 1-RDMs, i.e. the number of orbitals/sites
     */
    size_t dimension() const { return one_rdm.dimension(); }


    /**
     *  @return the difference between the alpha and beta 1-RDM
     */
    OneDM<Scalar> spinDensityRDM() const {
        return one_rdm_aa - one_rdm_bb;
    }
};


}  // namespace GQCP
