#include "euler_integrator.h"

DmsEuler::DmsEuler(PetscScalar dt, bool adpt, MPI_Comm comm) :
					DmsIntegrator::DmsIntegrator(dt, adpt, comm) {

	PetscFunctionBegin;
}

PetscErrorCode DmsEuler::integrate(const CVec& Coords, const CVec& RHS) {
	PetscFunctionBegin;

	for(int dim = 0; dim < Coords.size(); dim++) {
		ierr = VecAXPY(Coords[dim], Delta, RHS[dim]);
		DMS_CHKERRQ(ierr);
	}

	PetscFunctionReturn(ierr);
}
