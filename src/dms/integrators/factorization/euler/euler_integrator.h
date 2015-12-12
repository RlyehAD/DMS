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
 @Created: May 21, 2014

 * The DMS_Integraotr object simply advances the CG variables (corresponding to mesoscopic
 * or macroscopic classes) in time. These variables could be positions for the non-inertial
 * case, or both positions and velocities for the inertial one. Furthermore, the integrator
 * could be used to implement a (Trotter) factorization method, Langevin dynamics (that
 * results from multiscale perturbation), or etc. Independent of the method of choice, any
 * new integrator should publicly inherit from DMS_Integrator and override the method Integrate().
 *
 * The integrator class is completely side-effect free. Once instantiated, this class never
 * changes any of its properties but will only return the new CG state based on the supplied RHS
 * and timestep (Delta) of choice.
 *
 * TODO: Adaptive time stepping
*/

#ifndef DMS_EULER_H
#define DMS_EULER_H

#ifdef __cplusplus

#include "../../dms_integrator.h"

class DmsEuler : public DmsIntegrator {
public:
	explicit DmsEuler(PetscScalar step, bool adaptive = false, MPI_Comm = PETSC_COMM_SELF);
	virtual ~DmsEuler() {};
	virtual PetscErrorCode integrate(const CVec&, const CVec&);
};

#endif

#endif
