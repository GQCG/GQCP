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
#pragma once


namespace GQCP {


/**
 *  An elementary calculation that is regarded as one step in an algorithm.
 * 
 *  Derived classes should implement:
 *      - execute()
 * 
 *  @param _Environment             the type of the environment that this step can read from and write to
 */
template <typename _Environment>
class Step {
public:
    using Environment = _Environment;


public:

    /*
     *  DESTRUCTOR
     */

    virtual ~Step() = default;


    /*
     *  PUBLIC METHODS
     */

    /**
     *  Execute/perform this algorithm step.
     * 
     *  @param environment              the environment that this step can read from and write to
     */
    virtual void execute(Environment& environment) = 0;
};


}  // namespace GQCP
