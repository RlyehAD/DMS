#include "fvm.h"
#include "common.h"
#include "petscdt.h"

inline PetscScalar fvm::kernelFunction(vector<PetscScalar>& coords, const vector<PetscScalar> &gridPos, PetscScalar resol) {
	PetscScalar dotProd = .0;

	for(int d = 0; d < gridPos.size(); d++)
		dotProd += pow(gridPos[d] - coords[d], 2.0);

	dotProd = dotProd / pow(resol, 2.0);

	return exp(-dotProd);
}

vector<PetscScalar> fvm::coarseGrain(const Vec Coords, const Vec Coords_prev, DmsBase& Dbase, std::fstream& fpLog) {

	/***********************************************************************/
	/****** Computes the field variables of choice: MPI not YET DONE! *****/
	/*********************************************************************/

	PetscFunctionBegin;

	if(Dbase.getNcgAdj() <= 0)
		throw;

	auto fvList = Dbase.getFVlist(),
	     grid = Dbase.getGrid(),
	     it = fvList.begin(); // FVList contains atomic indices for each FV cell / node

	std::vector<PetscScalar> fieldVars(Dbase.getNcgAdj(), 0.f);
	vector<PetscScalar> gridPos(Dbase.getDim()), coords(Dbase.getDim());

	PetscScalar **nCoordsPtr = new PetscScalar[Dbase.getDim()];
	PetscScalar resol = Dbase.getResol();

	for(auto dim = 0; dim < Dbase.getDim(); dim++) {
  		ierr = VecGetArray(Coords()[dim], &(*nCoordsPtr+dim));
  		CHKERRQ(ierr);
  	}

	std::vector<PetscInt> indices(Dbase.getNcgAdj());
	for(auto i = 0; i < Dbase.getNcgAdj(); i++)
		indices[i] = i;

	for(auto i = 0; i < Dbase.getNcgAdj(); i++) {

		vector<PetscScalar> gridPos = grid[i];
		
		for(auto aa_it = it->begin(); aa_it != it->end(); aa_it++) {

			for(auto dim = 0; dim < Dbase.getDim(); dim++)
				coords[dim] = nCoords[dim][*aa_it];

			fieldVars[i] += fvm::kernelFunction(coords, grid[i], resol); // should we include the mass here? mass[*aa_it]
		}

		it++;
	}

	for(auto i = 0; i < Dbase.getNcgAdj(); i++)
		ierr = VecSetValues(fieldVars, Dbase.getNcgAdj(), indices.data(), fieldVars.data(), INSERT_VALUES);

	for(auto dim = 0; dim < Dbase.getDim(); dim++) {
  		ierr = VecRestoreArray(nCoords[dim], &(*nCoordsPtr+dim));
  		CHKERRQ(ierr);
  	}

  	delete[] nCoordsPtr;

	PetscFunctionReturn(ierr);
}
