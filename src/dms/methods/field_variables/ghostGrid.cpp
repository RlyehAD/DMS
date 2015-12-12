#include "fvm.h"

/* This file contains the functions needed in constructing, updating, and maintaining 
*  the ghost grid and all of its associated data structures. The ghiost grid is just like
*  the usual orthogonal grid except it is reduced in size to take into account only cells 
*  that slowly vary in time. Thus, the name "ghost".
*/

PetscErrorCode fvm::constructGridID(DmsBase& Dbase) {

  /******************************************************************************************/
  /****** The idea behind _GridID is similar to that used in finite difference methods *****/
  /****** Each grid cell is assigned an ID. This way, boundary conditions and adaptive ****/
  /****** Discretization become easy to implement and avoid an indexing nightmare.   *****/
  /**************************************************************************************/

  PetscFunctionBegin;
  PetscErrorCode ierr = 0;
  PetscInt count = 0;

  for(auto i = 0; i < Dbase.getNcgX() + Dbase.getNeighX(); i++)
    for(auto j = 0; j < Dbase.getNcgY() + Dbase.getNeighY(); j++)
      for(auto k = 0; k < Dbase.getNcgZ() + Dbase.getNeighZ(); k++)

        if(i == 0 | j == 0 | k == 0 | i == (Dbase.getNcgX() + 1)
           | j == (Dbase.getNcgY() + 1) | k == (Dbase.getNcgZ() + 1))
            gridID[i][j][k] = boundary;
        else
            gridID[i][j][k] = count++;

  PetscFunctionReturn(ierr);
}

PetscErrorCode fvm::constructGrid(const PetscInt NumNodes, DmsBase& Dbase) {

  /*************************************************************/
  /****** Discretize in space and construct a grid tensor *****/
  /***********************************************************/

  PetscFunctionBegin;
  PetscErrorCode ierr = 0;

  auto grid = Dbase.getGrid(),
       gridID = Dbase.getGridID(),
       box = Dbase.getBox(Dbase.getCoords());

  std::vector<PetscScalar> boxMin(Dbase.getDim());

  for(auto dim = 0; dim < Dbase.getDim(); dim++)
      VecMin(Dbsase.getCoords()[dim], NULL, &boxMin[dim]);

  if(grid.size() > 0)
	    grid.clear();

  grid.resize(Dbase.getNcg());
  // This is a 3 * NumNodes_x * NumNodes_y * NumNodes_z tensor

  for(auto i = 0; i < grid.size(); i++)
	  grid[i].resize(Dbase.getDim());


  PetscInt dimX = 0, dimY = 1, dimZ = 2;
  PetscScalar hx, hy, hz;

  for(auto i = 0; i < Dbase.getNcgX(); i++)
	for(auto j = 0; j < Dbase.getNcgY(); j++)
		for(auto k = 0; k < Dbase.getNcgZ(); k++) {

			hx = PetscScalar(i) / (Dbase.getNcgX() - 1.0) * box[dimX] + boxMin[dimX];
			hy = PetscScalar(j) / (Dbase.getNcgY() - 1.0) * box[dimY] + boxMin[dimY];
			hz = PetscScalar(k) / (Dbase.getNcgZ() - 1.0) * box[dimZ] + boxMin[dimZ];

			if(gridID[i][j][k] != boundary && gridID[i][j][k] != garbage) {
				grid[gridID[i][j][k]][dim_x] = hx;
				grid[gridID[i][j][k]][dim_y] = hy;
				grid[gridID[i][j][k]][dim_z] = hz;
			}
		}

  PetscFunctionReturn(ierr);
}
