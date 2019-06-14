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
#ifndef GQCP_NEWTONORBITALOPTIMIZATIONOPTIONS_HPP
#define GQCP_NEWTONORBITALOPTIMIZATIONOPTIONS_HPP


#include "OrbitalOptimization/OrbitalOptimizationOptions.hpp"

#include "math/optimization/BaseHessianModifier.hpp"
#include "math/optimization/UnalteringHessianModifier.hpp"

#include <utility>


namespace GQCP {


class NewtonOrbitalOptimizationOptions : public OrbitalOptimizationOptions {
private:
    std::shared_ptr<BaseHessianModifier> hessian_modifier;


public:
    // CONSTRUCTORS

    /**
     *  @param hessian_modifier         the modifier functor that should be used when an indefinite Hessian is encountered
     */
    NewtonOrbitalOptimizationOptions(std::shared_ptr<BaseHessianModifier> hessian_modifier = std::make_shared<UnalteringHessianModifier>(), const double convergence_threshold = 1.0e-08, const double maximum_number_of_iterations = 128);
};


}  // namespace GQCP



#endif  // GQCP_NEWTONORBITALOPTIMIZATIONOPTIONS_HPP
