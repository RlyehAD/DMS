#include "fvm.h"
#include "common.h"

/* Data structures for FVM
**
*/

PetscErrorCode fvm::ConstructFVList(DmsBase& Dbase) {

  /*****************************************************************************/
  /****** Setup node indices list for efficient computation of the Kernel *****/
  /****** Here the GridID is updated to take into account the transforma- ****/
  /****** tion of the Grid to a ghost Grid that excludes certain nodes   ****/
  /*************************************************************************/
  PetscErrorCode ierr = 0;
  auto fvList = Dbase.getFVlist();
  auto gridID = Dbase.getGridID();

  fvList.clear();
  fvList.resize(Dbase.getNcg());

  auto box = Dbase.getBox(Dbase.Microscopic.getCoords());
  CVec nCoords(Dbase.getDim());

  for(auto dim = 0; dim < Dbase.getDim(); dim++) {

  		ierr = VecCreateMPI(Dbase.getComm(), PETSC_DECIDE, Dbase.getNatoms(), nCoords.data() + dim);
  		CHKERRQ(ierr);

  		ierr = VecCopy(Dbase.Microscopic.getCoords()[dim], nCoords[dim]);
  		CHKERRQ(ierr);

		ierr = VecAssemblyBegin(nCoords[dim]);
		CHKERRQ(ierr);
		ierr = VecAssemblyEnd(nCoords[dim]);
		CHKERRQ(ierr);

  		ierr = VecScale(nCoords[dim], box[dim]);
  		CHKERRQ(ierr);
		ierr = VecAssemblyBegin(nCoords[dim]);
		CHKERRQ(ierr);
		ierr = VecAssemblyEnd(nCoords[dim]);
		CHKERRQ(ierr);
  }

  PetscScalar **nCoordsPtr = new PetscScalar[Dbase.getDim()];

  for(auto dim = 0; dim < Dbase.getDim(); dim++) {
  		ierr = VecGetArray(nCoords[dim], &(*nCoordsPtr+dim));
  		CHKERRQ(ierr);
  }

  PetscInt indX, indY, indZ;

  try {
	  for(auto i = 0; i < Dbase.getNatoms(); i++) {
		  PetscInt count = 0;

		  vector<PetscInt> atomic_indices(count);
		  count = 0;

	  	  indX = static_cast<PetscInt>(nCoordsPtr[0] / Box[0]);
	  	  indY = static_cast<PetscInt>(nCoordsPtr[1] / Box[1]);
	  	  indZ = static_cast<PetscInt>(nCoordsPtr[2] / Box[2]);

		  fvList[gridID[indX][indY][indZ]].push_back(i);
	  }
  }
  catch(int e)
	std::cout << "An exception occurred. Exception Nr. " << e << std::endl 
	<< "You most probably ran out of memory." << std::endl;

  for(auto dim = 0; dim < Dbase.getDim(); dim++) {
  		ierr = VecRestoreArray(nCoords[dim], &(*nCoordsPtr+dim));
  		CHKERRQ(ierr);
  }

  delete[] nCoordsPtr;

  PetscInt count = 0;
  auto it = fvList.begin();

  for(auto i = 0; i < Dbase.getNcgX() + Dbase.getNeighX(); i++)
	  for(auto j = 0; j < Dbase.getNcgY() + Dbase.getNeighY(); j++)
		  for(auto k = 0; k < Dbase.getNcgZ() + Dbase.getNeighZ(); k++)
			  if(gridID[i][j][k] >= 0) {
				  if(it->size() < Dbase.getThresh()) {
					  gridID[i][j][k] = garbage;
					  it = fvList.erase(it); // here the linked list destructor should be invoked
					  count++;
				  }
				  else
					  it++;
			  }

  Dbase.setNcgAdj(Dbase.getNcg() - count);

  if(count > 0) {

	  std::cout << "Based on the threshold supplied: " << Dbase.getThresh() << ", the number of nodes has" <<
		  " changed from " << Dbase.getNcg() << " to " << Dbase.getNcg() - count << std::endl;
	  // Remap GridID in case count is NOT zero
	  count = 0;

	  for(auto i = 0; i < Dbase.getNcgX() + Dbase.getNeighX(); i++)
		  for(auto j = 0; j < Dbase.getNcgY() + Dbase.getNeighY(); j++)
			  for(auto k = 0; k < Dbase.getNcgZ() + Dbase.getNeighZ(); k++)
				  if(gridID[i][j][k] != garbage && gridID[i][j][k] != boundary)
					  gridID[i][j][k] = count++;
  }

  PetscFunctionReturn(ierr);
}

PetscErrorCode fvm::deleteJacobians(DmsBase& Dbase) {
	PetscFunctionBegin;
	PetscErrorCode ierr;

	auto jacobKernel = getMesoMicroVec();

	for(auto dim = 0; dim < Dbase.getDim(); dim++) {
		ierr = MatDestroy(&(jacobKernel[dim]));
		CHKERRQ(ierr);
	}

	jacobKernel.clear(); // Very important for push pack to work when updating grid!!

	ierr = MatDestroy(jacobian);
	CHKERRQ(ierr);

	PetscFunctionReturn(ierr);
}

PetscErrorCode fvm::setupJacobians(DmsBase& Dbase) {

	/******************************************************************/
	/****** Allocates memory for the Jacobians: runs in MPI mode *****/
	/****************************************************************/

	/* This function allocates memory for the kernel and full Jacobian matrices.
	It does not do any computations. Therefore, when updating the grid or
	discretizing, it is absolutely *necessary* to call deleteJacobians()
	in order to avoid memory lacks.
	*/

	PetscFunctionBegin;
	PetscErrorCode ierr;

	auto fvList = Dbase.getFVlist(),
	     jacobKernel = getMesoMicroVec(),
	     jacobian = Dbase.getMesoMicro(),
	     gridID = getGridID();

	PetscInt *nnz = new PetscInt[Dbase.getNcgAdj()];
	PetscInt count = 0;

	for(auto it = fvList.begin(); it != fvList.end(); it++)
		nnz[count++] = it->size();

	for(auto d = 0; d < Dbase.getDim(); d++) {
		Mat matTmp;
		ierr = MatCreate(Dbase.getComm(), &matTmp);
		MatSetSizes(matTmp, PETSC_DECIDE, PETSC_DECIDE, Dbase.getNcgAdj(), Dbase.getNatoms()); // delete and reallocate when updating grid
		MatSetType(matTmp, MATMPIAIJ);
		MatMPIAIJSetPreallocation(matTmp, 0, nnz, 0, nnz);
		jacobKernel.push_back(matTmp);
		//MatView(FieldVar::JacobKernel[d], PETSC_VIEWER_STDOUT_SELF);
	}

	delete[] nnz;

	nnz = new PetscInt[Dbase.getNcgAdj()];

	for(auto i = 0; i < Dbase.getNcgAdj(); i++)
		nnz[i] = 0;

	// 1 ... (NumNodes + Grid_Neighbors - 1) -> exclude boundary points
	for(auto i = 1; i < Dbase.getNcgX() + Dbase.getNeighX() - 1; i++)
		for(auto j = 1; j < Dbase.getNcgY() + Dbase.getNeighY() - 1; j++)
			for(auto k = 1; k < Dbase.getNcgZ() + Dbase.getNeighZ() - 1; k++)
				if(gridID[i][j][k] >= 0 ) {

					PetscInt index_x_min = -1,
						 index_y_min = -1,
						 index_z_min = -1, // extended gird searching must be taken into account!!! [FUTURE]
						 index_x_max = +1,
						 index_y_max = +1,
						 index_z_max = +1;

					// neighbor searching along the x, y, and z- axes
					for(auto index_x = index_x_min; index_x <= index_x_max; index_x++)
						for(auto index_y = index_y_min; index_y <= index_y_max; index_y++)
							for(auto index_z = index_z_min; index_z <= index_z_max; index_z++)
								if(gridID[index_x+i][index_y+j][index_z+k] >= 0)
									nnz[gridID[i][j][k]] += 1;
				}

	ierr = MatCreate(Dbase.getComm(), jacobian));
	CHKERRQ(ierr);

	ierr = MatSetSizes(jacobian, PETSC_DECIDE, PETSC_DECIDE, Dbase.getNcgAdj(), Dbase.getNcgAdj());
	CHKERRQ(ierr);

	ierr = MatSetType(jacobian, MATMPIAIJ);
	CHKERRQ(ierr);

	ierr = MatMPIAIJSetPreallocation(jacobian, 0, nnz, 0, nnz);
	CHKERRQ(ierr);

	PetscFunctionReturn(ierr);
}

