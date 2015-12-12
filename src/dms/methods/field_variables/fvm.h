/*
 * This file is part of the DMS molecular simulation package.
 *
 * Copyright (c) 2014, by the DMS development team:,
 * Andrew Abi Mansour
 *
 * DMs is free software; you can redistribute it and/or
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
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund DMS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.dms.org

 @Author: Andrew Abi Mansour
 @Created: March 27, 2014

 *This file should include ONLY the functions to be interfaced between DMS (CG phase) and GROMACS (MD phase)
 *It should contain all the necessary functions related to the space warping method.
 *In order to interface classes with C, a void CSpaceWarping pointer is created to handle
 *all the possible methods to be called in md.c

 *The SpaceWarping class changes the (micro) state from MD. Essentially, it advances the atomic positions
 *in time via a CG-guided evolution.
*/

#ifndef DMS_FVM_H
#define DMS_FVM_H

class DmsBase;
enum cellTYpe {boundary = -1, garbage = -2};

namespace fvm {
	PetscErrorCode initialize(DmsBase& Dbase);	
	PetscErrorCode updateRef(DmsBase& Dbase);
	PetscErrorCode coarseGrain(const Vec Coords, const Vec pCoords, const PetscInt dim, DmsBase& Dbase);
	PetscErrorCode fineGrain(const Vec Coords, const Vec Coords_prev, const PetscInt dim, DmsBase& Dbase);
	PetscErrorCode constructVelocities(DmsBase& Dbase);
	PetscErrorCode constructCoords(const CVec& Vars, const CVec& pVars, DmsBase& Dbase);
}

#endif
