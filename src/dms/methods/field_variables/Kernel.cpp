#include "Interface.h"

PetscScalar FieldVar::KernelFunction(const PetscScalar* const Coords, const vector<PetscScalar> &GridPos) {
	PetscScalar dotprod = .0;

	for(int d = 0; d < FieldVar::_Dim; d++)
		dotprod += pow((GridPos[d] - Coords[d]) / FieldVar::Resol, 2.0);

	return exp(-dotprod);
}

vector<PetscScalar> FieldVar::ComputeFV(const PetscScalar* const Coords) {

	/***********************************************************************/
	/****** Computes the field variables of choice: MPI not YET DONE! *****/
	/*********************************************************************/

	PetscFunctionBegin;
	PetscInt count = 0;

	for(auto i = 0; i < FieldVar::NumNodes_x + FieldVar::Grid_Neighbors_x; i++)
		for(auto j = 0; j < FieldVar::NumNodes_y + FieldVar::Grid_Neighbors_y; j++)
			for(auto k = 0; k < FieldVar::NumNodes_z + FieldVar::Grid_Neighbors_z; k++) {

				PetscInt grid_id = FieldVar::GridID[i][j][k];
				if(grid_id >= 0) {
					FieldVar::FieldVars[count] = 0.f;
					vector<PetscScalar> GridPos = FieldVar::Grid[grid_id];

					for(auto it = FieldVar::AtomicList[grid_id].begin(); it != FieldVar::AtomicList[grid_id].end(); it++)
						FieldVar::FieldVars[count] += FieldVar::KernelFunction(Coords+(FieldVar::Dim)*(*it), GridPos);

					count++;
				}
			}

	 return FieldVar::FieldVars;
}

void FieldVar::Py_ComputeCG_Pos(double *COORDS_IN, int NATOMS, int DIM, double *CG_OUT, int NUMCG) {
	PetscFunctionBegin;
	vector<PetscScalar> tmp = FieldVar::_ComputeFV(COORDS_IN);

	for(auto i = 0; i < tmp.size(); i++)
		CG_OUT[i] = tmp[i];

}

vector<PetscScalar> FieldVar::ComputeFV_Vel(const PetscScalar* const Coords, const PetscScalar* const Vel) {

	/*************************************************************************/
	/****** Computes the velocity of field variables: MPI not YET DONE! *****/
	/***********************************************************************/

	PetscFunctionBegin;
	PetscInt count = 0;

	for(auto i = 0; i < FieldVar::NumNodes_x + FieldVar::Grid_Neighbors_x; i++)
		for(auto j = 0; j < FieldVar::NumNodes_y + FieldVar::Grid_Neighbors_y; j++)
			for(auto k = 0; k < FieldVar::NumNodes_z + FieldVar::Grid_Neighbors_z; k++) {

				PetscInt grid_id = FieldVar::GridID[i][j][k];
				if(grid_id >= 0) {
					FieldVar::FieldVars_Vel[count] = 0.f;
					vector<PetscScalar> GridPos = FieldVar::Grid[grid_id];

					for(auto it = FieldVar::AtomicList[grid_id].begin(); it != FieldVar::AtomicList[grid_id].end(); it++) {

						const PetscScalar *const AtomicCoords = Coords + (FieldVar::Dim)*(*it),
								  *const AtomicVel = Vel + (FieldVar::Dim)*(*it),
								   Kernel = FieldVar::KernelFunction(AtomicCoords, GridPos);

						for(auto d = 0; d < FieldVar::Dim; d++)
							FieldVar::FieldVars_Vel[count] += 2.0 * (GridPos[d] - AtomicCoords[d]) / pow(FieldVar::Resol,2.0) * AtomicVel[d] * Kernel;
					}
					count++;
				}
			}

	 return FieldVar::FieldVars_Vel;
}

void FieldVar::Py_ComputeCG_Vel(double *COORDS_IN, int NATOMS1, int DIM1, double *VEL_IN, int NATOMS2, int DIM2, double *CG_OUT, int NUMCG) {
	PetscFunctionBegin;
	vector<PetscScalar> tmp = FieldVar::ComputeFV_Vel(COORDS_IN, VEL_IN);

	for(auto i = 0; i < tmp.size(); i++)
		CG_OUT[i] = tmp[i];
}

vector<PetscScalar> FieldVar::ComputeFV_For(const PetscScalar* const Coords, const PetscScalar* const Vel, const PetscScalar* const For) {

	/**********************************************************************/
	/****** Computes the force of field variables: MPI not YET DONE! *****/
	/********************************************************************/

	PetscFunctionBegin;
	PetscInt count = 0;

	for(auto i = 0; i < FieldVar::NumNodes_x + FieldVar::Grid_Neighbors_x; i++)
		for(auto j = 0; j < FieldVar::NumNodes_y + FieldVar::Grid_Neighbors_y; j++)
			for(auto k = 0; k < FieldVar::NumNodes_z + FieldVar::Grid_Neighbors_z; k++) {

				PetscInt grid_id = FieldVar::GridID[i][j][k];
				if(grid_id >= 0) {
					FieldVar::_FieldVars_For[count] = 0.f;
					vector<PetscScalar> GridPos = FieldVar::Grid[grid_id];

					for(auto it = FieldVar::AtomicList[grid_id].begin(); it != FieldVar::AtomicList[grid_id].end(); it++) {

						const PetscScalar *const AtomicCoords = Coords+(FieldVar::_Dim)*(*it),
										  * const AtomicVel = Vel+(FieldVar::_Dim)*(*it),
										  * const AtomicFor = For+(FieldVar::_Dim)*(*it),
										  Kernel = FieldVar::_KernelFunction(AtomicCoords, GridPos);

						for(auto d = 0; d < FieldVar::_Dim; d++)
							FieldVar::_FieldVars_For[count] += 2.0 * ( (GridPos[d] - AtomicCoords[d]) / pow(FieldVar::_Resol,2.0) * AtomicFor[d]
							          - AtomicVel[d] / pow(FieldVar::_Resol,2.0) * AtomicVel[d]
							          + 2.0 * AtomicVel[d] * pow((GridPos[d] - AtomicCoords[d]) / pow(FieldVar::_Resol,2.0), 2.0) * AtomicVel[d] )
							          * Kernel;
					}
					count++;
				}
			}

	 return FieldVar::_FieldVars_For;
}

void FieldVar::Py_ComputeCG_For(double *COORDS_IN, int NATOMS1, int DIM1, double *VEL_IN, int NATOMS2, int DIM2, double *FOR_IN, int NATOMS3, int DIM3, double *CG_OUT, int NUMCG) {
	PetscFunctionBegin;
	vector<PetscScalar> tmp = FieldVar::_ComputeFV_For(COORDS_IN, VEL_IN, FOR_IN);

	for(auto i = 0; i < tmp.size(); i++)
		CG_OUT[i] = tmp[i];
}
