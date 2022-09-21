#include "dms_integrator.h"

DmsIntegrator::DmsIntegrator(PetscScalar dt, bool adpt, MPI_Comm COMM) : Delta(dt),
								   adaptive(adpt), comm(COMM) {

	PetscFunctionBegin;
	PetscFunctionReturn(0);
}
