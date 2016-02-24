 /*
 * This file is part of the DMS molecular simulation package.
 *
 * @Author: Andrew Abi Mansour
 * @Created: March 24, 2014
 *
 * DMS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * DMS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with DMS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to DMS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official DMS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.
 *
 * To help us fund DMS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.
 *
 * About this file:

 *This file should include ONLY the functions to be interfaced between DMS (CG phase) and GROMACS (MD phase)
 *It should include all header files specific to the coarse graining method (space warping, field variables, etc.)
 *Each CG method method is represented by an object that *must* include the following methods:

 -CoarseGrain :  reduce the dimensionality by constructing CGs, etc.
 -Integrate   :  advance the CG state in time
 -FineGrain   :  recover the atomistic state

 *The three functions will be called from one interface function (CG_step) that should available as a public method.
 *The class must be instantiated before the multiscale simulation begins and must persist throughout the entire period.
 *The atomistic and CG states must also persist throughout the simulation. Both are updated via the class methods.
 *Therefore, all three methods need not return anything.

 */

#ifndef DMS_H
#define DMS_H

#include "types/enums.h" // must include this first for state.h
#include "types/state.h"
#include "gmx_fatal.h"

#ifdef __cplusplus

#include <iostream>
#include <vector>
#include <petsc.h>
#include <petscmat.h>
#include <petscsys.h>
#include <memory>

extern "C" {
#endif

#include "methods/space_warping/interface_swm.h"

#ifdef __cplusplus
}
#endif

#endif
