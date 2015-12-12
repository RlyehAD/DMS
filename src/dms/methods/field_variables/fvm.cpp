#include "fvm.h"

/* This file contains some of the functions pertaining to the field variables method (FVM).
** TODO: I need to to put more info here.
** 
**
**/

PetscErrorCode fvm::initialize(DmsBase& Dbase)

	PetscFunctionBegin;

	// GridID must remain constant throughout the simulation
	// GridID keeps track of grid points that should be included
	auto gridID = Dbase.getGridID();

	gridID.resize(Dbase.getNcgX() + Dbase.getNeighX()); // why do we include neighbor cells here?!!!!!!!!!!!
	for(auto i = 0; i < gridID.size(); i++) {
		gridID[i].resize(Dbase.getNcgY() + Dbase.getNeighY());
		for(auto j = 0; j < gridID[i].size(); j++)
			gridID[i][j].resize(Dbase.getNcgZ()+ Dbase.getNeighZ());
	}

	// structured grid IDs all set
	ierr = fvm::constructGridID(Dbase);
	CHKERRQ(ierr);

	ierr = fvm::constructGrid(Dbase.getNcg(), Dbase);
	CHKERRQ(ierr);

  	// AdjNumNodes setup here!
  	ierr = fvm::constructFVList(Dbase); 
  	CHKERRQ(ierr);

  	// reconstruct the ghost grid
  	ierr = fvm::constructGrid(Dbase.getNcgAdj(), Dbase); 
  	CHKERRQ(ierr);

  	ierr = fvm::setupJacobians();
  	CHKERRQ(ierr);

  	PetscFunctionReturn(ierr);
}

PetscErrorCode SpaceWarping::updateRef() {

	if(!(TimeStep % FreqUpdate)) {
        	// Microscopic->Ref_Coords is updated here only
                std::cout << "Updating reference structure ..." << std::endl;
                ierr = computeKernelJacob(Microscopic->Get_Coords());
		// must also update cons Jacob
                CHKERRQ(ierr);
	}

	PetscFunctionReturn(ierr);
}
