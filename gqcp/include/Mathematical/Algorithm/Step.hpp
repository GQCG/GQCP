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


#include <string>


namespace GQCP {


/**
 *  An elementary calculation that is regarded as one step in an algorithm.
 * 
 *  Derived classes should implement:
 *      - description()
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
     *  VIRTUAL PUBLIC METHODS
     */

    /**
     *  @return a textual description of this algorithmic step
     */
    virtual std::string description() const = 0;

    /**
     *  Execute/perform this algorithm step.
     * 
     *  @param environment              the environment that this step can read from and write to
     */
    virtual void execute(Environment& environment) = 0;
};


}  // namespace GQCP
