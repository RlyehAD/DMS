#include "pade_integrator.h"

DmsPade::DmsPade(PetscScalar dt, PetscScalar dtMicro, bool Adpt, MPI_Comm comm) :
  DmsIntegrator::DmsIntegrator(dt, Adpt, comm), delta(dtMicro) {

	PetscFunctionBegin;

}

PetscErrorCode DmsPade::computeInvMat(PetscScalar Coords[]) {
  /* This function computes the inverse matrix of the PAs for computing the coefficients
     a1, a1/2, and b1. See ref.
  */
        PetscFunctionBegin;
/*
	PetscScalar coefsRaw[nHist*nHist] = { 
	  delta, pow(delta, 0.5), - Coords[0] * delta, 
	  Delta, pow(Delta, 0.5), - Coords[1] * Delta, 
	  Delta + delta, pow(Delta + delta, 0.5), - Coords[2] * (Delta + delta)
	};

	ierr = MatCreateSeqDense(comm, nHist, nHist, coefsRaw, &invMat);
	CHKERRQ(ierr);

	I commented this block because gcc/4.9.2 was complaining about uninitialized coefsRaw. WTF??
*/

	PetscFunctionReturn(ierr);
}

PetscErrorCode DmsPade::integrate(const std::vector< CVec >& Coords, const CVec& RHS) {
	PetscFunctionBegin;

	std::cout << "Calling DmsPade ..." << std::endl;

        for(int dim = 0; dim < Coords[1].size(); dim++) {
                ierr = VecAXPY(Coords[1][dim], delta, RHS[dim]);
                DMS_CHKERRQ(ierr);
        }

	/*
	for(auto dim = 0; dim < Coords[0].size(); dim++)
		for(auto order = 0; order < Coords[0][0].size(); dim++) { 
		
			PetscScalar CoordsTmp[] = {Coords[0][dim][order], Coords[0][dim][order], Coords[0][dim][order]};
			ierr = VecAXPY(CoordsTmp, Delta, RHS[dim]);
			CHKERRQ(ierr);
		}
	}
	*/

	PetscFunctionReturn(ierr);
}
