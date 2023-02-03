#include "common.h"
#include "swm.h"

/* The SWM class is created in serial on all processors, so the communicator is naturally
 * MPI_COMM_SELF. The reason is SWM uses a small number of CG variables to capture the overall
 * global changes of a given macromolecule.
 *
 * TODO: Support subsystem decomposition
 */

PetscErrorCode swm::initialize(DmsBase& Dbase) {

	PetscFunctionBegin;
	PetscErrorCode ierr;

	std::vector< std::vector<PetscInt> > indices = computeIndices(Dbase.getOrder());

	Dbase.setNcg(indices.size());

	ierr = setupBasis(Dbase.getLocalNatoms(), Dbase.getNcg(), Dbase);
	CHKERRQ(ierr);

	ierr = constructBasis(Dbase.Microscopic->Get_RefCoords(), Dbase, Dbase.fpLog);
	CHKERRQ(ierr);

	PetscFunctionReturn(ierr); // this is redundant, indices are being computed elsewhere
}

std::vector< std::vector<PetscInt> > swm::computeIndices(int kmax) {
	PetscFunctionBegin;
	std::vector< std::vector<PetscInt> > indices_tmp;

	for (int n = 0; n < kmax + 1; n++)
	       	for (int i = 0; i < n + 1; i++)
	          	for (int j = 0; j < n + 1-  i; j++) {
	           		PetscInt tmp[] = {n - i -j, j, i};
	            		std::vector<PetscInt> tmpv(tmp, tmp + sizeof(tmp) / sizeof(tmp[0]));

				//if(!(n == 0 && i == 0 && j == 0))
	            			indices_tmp.push_back(tmpv);
	            	}

	PetscFunctionReturn(indices_tmp);
}

PetscErrorCode swm::setupBasis(const PetscInt nAtoms, const PetscInt nCG, DmsBase& Dbase) {
	PetscFunctionBegin;

	PetscErrorCode ierr;
	Mat *mesoMicroMap = Dbase.getMesoMicro(),
	    *microMesoMap = Dbase.getMicroMeso(),
	    *kernel = Dbase.getKernel();
 
	ierr = MatCreateSeqDense(Dbase.getComm(), nAtoms, nCG, NULL, mesoMicroMap);
	CHKERRQ(ierr);

	ierr = MatCreateSeqDense(Dbase.getComm(), nAtoms, nCG, NULL, kernel);
	CHKERRQ(ierr);

	ierr = MatCreateSeqDense(Dbase.getComm(), nCG, nCG, NULL, microMesoMap);
	CHKERRQ(ierr);

	PetscFunctionReturn(ierr);
}

PetscErrorCode swm::updateRef(DmsBase& Dbase) {

	PetscErrorCode ierr;

	if(!(Dbase.getTimeStep() % Dbase.getFreqUpdate())) {

        	// Microscopic->Ref_Coords is updated here only
                ierr = constructBasis(Dbase.Microscopic->Get_RefCoords(), Dbase, Dbase.fpLog);
                CHKERRQ(ierr);
		

		for(int dim = 0; dim < Dbase.Mesoscopic->Get_Dim(); dim++) {
                        ierr = VecCopy(Dbase.Mesoscopic->Get_RefCoords()[dim], Dbase.Mesoscopic->Get_Coords()[dim]);
                        CHKERRQ(ierr);
                }

	}

	PetscFunctionReturn(ierr);
}

PetscErrorCode swm::constructCoords(CVec Vars, CVec pVars, DmsBase& Dbase, std::fstream& fp) {
	/* This function constructs SWM coords, velocities, and forces, depending on the arguments
 	* supplied. The three CG-construction function pointers in Mesoscopic should all point to
	 * this function.
	 */
	PetscFunctionBegin;
	PetscErrorCode ierr;

	ierr = (Dbase.Mesoscopic->mapping)(Vars, pVars, Dbase, fp);
	CHKERRQ(ierr);

	PetscFunctionReturn(ierr);
}

PetscErrorCode swm::constructVelocities(DmsBase& Dbase) {
	PetscFunctionBegin;
	PetscErrorCode ierr;

	for(int dim = 0; dim < Dbase.Mesoscopic->Get_Dim(); dim++) {

		ierr = VecCopy(Dbase.Mesoscopic->Get_Coords()[dim],
			       Dbase.Mesoscopic->Get_Velocities()[dim]);
		CHKERRQ(ierr);

		ierr = VecAXPY(Dbase.Mesoscopic->Get_Velocities()[dim], -1.0,
			       Dbase.Mesoscopic->Get_pCoords()[dim]);
		CHKERRQ(ierr);

		ierr = VecScale( Dbase.Mesoscopic->Get_Velocities()[dim], 1.0 / Dbase.Microscopic->getLength() );

		CHKERRQ(ierr);
	}

	PetscFunctionReturn(ierr);
}
