// This file is part of mpc.

// mpc is free software: you can redistribute it and/or
// modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// mpc is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with mpc.  If not, see
// <http://www.gnu.org/licenses/>.

#pragma once

// stl
#include <iostream>

#ifndef _DEBUG
#define CONSTRAINT_DELETION_WARN(warn, format, ...)
#else
#define CONSTRAINT_DELETION_WARN(warn, format, ...) \
    if ((warn)) {                                   \
        fprintf(stderr, format, __VA_ARGS__);       \
    }
#endif