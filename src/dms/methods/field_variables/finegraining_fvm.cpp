#include "fvm.h"
#include "common.h"

PetscErrorCode fvm::fineGrain(const Vec Coords, const Vec Coords_prev, DmsBase& Dbase, std::fstream& fpLog)
{
	/* Finegraining method: r = r_prev + U * (\phi - \phi_prev) + \sigma
	 * \sigma is just the residuals vector, which we can compute iteratively but
	 * for simplicity is taken to be zero here.
	 */
	PetscFunctionBegin;
	// Begin Newton && steepest descent method

	fpLog << Dbase.getTime() << ":INFO:Attempting to recover the microstate iteratively" << std::endl;
 
	Vec atomicDisp;
	VecCreateMPI(Dbase.getComm(), PETSC_DECIDE, Dbase.getNatoms(), &atomicDisp);

	PetscInt iters = 0;
	PetscScalar consError;

	while(atomicError >= Dbase.getTol()) {

		PetscScalar scaling = max(Dbase.getScale() / (iters*0.1 + 1.0), 0.1);
		consError = fvm::computeLagrangeMulti(Dbase.Microscopic.getCoords(), multipliers, Dbase.Mesoscopic.getCoords()[0], scaling, iters, fpLog, Dbase.getAssemFreq());

		PetscScalar multiNorm;
		atomicError = .0;

		for(auto dim = 0; dim < Dbase.getDim(); dim++) {
			Vec tmp = fvm::advancedAtomistic(Dbase.Microscopic.getCoords()[dim], multipliers, dim);

			VecCopy(tmp, atomicDisp);
			VecAXPBY(atomicDisp, -1.0, 1.0, Dbase.Microscopic.getCoords()[dim]); // atom_disp -= Coords_petsc

			PetscScalar atomicTmp;
			VecNorm(atomicDisp, NORM_INFINITY, &atomicTmp);
			VecCopy(tmp, Dbase.getCoords()[dim]);

			VecDestroy(&tmp);
			atomicError = max(atomicTmp, atomicError);
		}

		iters++;
		fpLog << Dbase.getTime() << ":INFO:Newton max atomic displacement: " << atomic_error << std::endl;
		fpLog << Dbase.getTime() <<":INFO:Constraint error = " << cons_error << endl;
	}

	PetscFunctionReturn(ierr);
}

PetscScalar fvm::computeLagrangeMulti(CVec coords, Vec multipliers, Vec fieldVars, PetscScalar scaling, PetscInt new_iters, std::fstream fpLog, bool assemble) {
	PetscFunctionBegin;

	PetscErrorCode ierr;
	auto jacobian = Dbase.getJAcobian();

	Vec cons = constraints(coords);
	fpLog << FieldVar::getTime() << ":INFO:Computing Lagrange multipliers" << endl;

	VecAXPBY(cons, 1.0, -1.0, fieldVars); // Cons = FV - Cons

	PetscScalar cons_error;
	VecNorm(cons, NORM_INFINITY, &cons_error);

	if(assemble)
		if(new_iters < 1) {
			ierr = assembleJacobian(coords);
			CHKERRQ(ierr);

			ierr = MatShift(jacobian, scaling);
			CHKERRQ(ierr);
		}

	//MatView(FieldVar::Jacobian, PETSC_VIEWER_STDOUT_WORLD);

	fpLog << Dbase.getTime() << ":INFO:Calling KSP solver " << std::endl;

	KSP ksp;
	PetscInt iters;

	ierr = KSPCreate(Dbase.getTime(), &ksp);
	CHKERRQ(ierr);

	ierr = KSPSetOperators(ksp, jacobian, jacobian, SAME_PRECONDITIONER);
	CHKERRQ(ierr);
	ierr = KSPSetType(ksp, KSPBCGS);
	CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp);

	ierr = KSPSolve(ksp, cons, multipliers);
	CHKERRQ(ierr);
	ierr = KSPGetIterationNumber(ksp, &iters);
	CHKERRQ(ierr);
	ierr = KSPDestroy(&ksp);
	CHKERRQ(ierr);

	ierr = VecDestroy(&cons);
	CHKERRQ(ierr);

	fpLog << Dbase.getTime() << ":INFO:KSP converged in " << iters << std::endl;

	PetscFunctionReturn(cons_error);
}

PetscErrorCode fvm::kernelJacobian(const Vec* const Coords, const PetscInt &dim) {

	/******************************************************************************/
	/******** Construct the transpose of the Jacobian of the kernel function *****/
	/************ The size of the Jacobian is Natoms x N_CG  ********************/
	/***************************************************************************/

	/* The kernel Jacobian is NOT initialized to zero every time it is updated
	*  because the nnz entries are INSERTED at the same location every time. 
	*/
	
	PetscFunctionBegin;

        FieldVar::fp << FieldVar::GetTime() << ":INFO:Computing Kernel Jacobian " << endl;

	PetscErrorCode ierr;
	PetscInt istart, iend;

	VecGetArray(Coords[0], &(FieldVar::Coords_Local_x));
	VecGetArray(Coords[1], &(FieldVar::Coords_Local_y));
	VecGetArray(Coords[2], &(FieldVar::Coords_Local_z));

	MatGetOwnershipRange(FieldVar::JacobKernel[dim], &istart, &iend);

	for(auto i = istart; i < iend; i++) {

		PetscScalar *Jk = new PetscScalar[FieldVar::FVList[i].size()];
		PetscInt *Indices = new PetscInt[FieldVar::FVList[i].size()];
		PetscInt count = 0;
				    
		for(auto atom = FieldVar::FVList[i].begin(); atom != FieldVar::FVList[i].end(); atom++) {
			vector<PetscScalar> GridPos = FieldVar::Grid[i];
			const PetscScalar r[] = {FieldVar::Coords_Local_x[*atom], FieldVar::Coords_Local_y[*atom], FieldVar::Coords_Local_z[*atom]};
			Jk[count] = (FieldVar::Grid[i][dim] - r[dim]) * (FieldVar::KernelFunction(r, FieldVar::Mass[*atom], GridPos)) * (2.0 / pow(FieldVar::Resol,2.0));
			Indices[count++] = (*atom);
		}

		ierr = MatSetValues(FieldVar::JacobKernel[dim], 1, &i, FieldVar::FVList[i].size(), Indices, Jk, INSERT_VALUES); CHKERRQ(ierr);
		delete[] Jk;
		delete[] Indices;
	}

	ierr = MatAssemblyBegin(FieldVar::JacobKernel[dim], MAT_FINAL_ASSEMBLY);
	ierr = MatAssemblyEnd(FieldVar::JacobKernel[dim], MAT_FINAL_ASSEMBLY);

	//MatView(FieldVar::JacobKernel[dim], PETSC_VIEWER_STDOUT_SELF);
	VecRestoreArray(Coords[0], &(FieldVar::Coords_Local_x));
	VecRestoreArray(Coords[1], &(FieldVar::Coords_Local_y));
	VecRestoreArray(Coords[2], &(FieldVar::Coords_Local_z));

	PetscFunctionReturn(ierr);
}

PetscErrorCode fvm::assembleJacobian(const Vec *const Coords) {
	PetscFunctionBegin;

        FieldVar::fp << FieldVar::GetTime() << ":INFO:Assembling Jacobian " << endl;

	PetscErrorCode ierr;
	auto node = (FieldVar::FVList).begin();

        VecGetArray(Coords[0], &(FieldVar::Coords_Local_x));
        VecGetArray(Coords[1], &(FieldVar::Coords_Local_y));
        VecGetArray(Coords[2], &(FieldVar::Coords_Local_z));

	for(auto i = 1; i < FieldVar::NumNodes_x + FieldVar::Grid_Neighbors_x - 1; i++)
		for(auto j = 1; j < FieldVar::NumNodes_y + FieldVar::Grid_Neighbors_y - 1; j++)
			for(auto k = 1; k < FieldVar::NumNodes_z + FieldVar::Grid_Neighbors_z - 1; k++) {

				PetscInt grid_index1 = (FieldVar::GridID)[i][j][k];

				if(grid_index1 >= 0) { // ghost cells and boundary cells not taken into account

					PetscInt index_x_min = -1,
						 index_y_min = -1,
						 index_z_min = -1, // extended gird searching must be taken into account!!! [FUTURE]
                	             		 index_x_max = +1,
                        	     		 index_y_max = +1,
                             			 index_z_max = +1;

					// Check for boundary cells
					if(FieldVar::GridID[i-1][j][k] < 0)
						index_x_min = 0;

					if(FieldVar::GridID[i+1][j][k] < 0)
						index_x_max = 0;

					if(FieldVar::GridID[i][j-1][k] < 0)
						index_y_min = 0;

					if(FieldVar::GridID[i][j+1][k] < 0)
						index_y_max = 0;

					if(FieldVar::GridID[i][j][k-1] < 0)
						index_z_min = 0;

					if(FieldVar::GridID[i][j][k+1] < 0)
						index_z_max = 0;

					vector<PetscScalar> GridPos = (FieldVar::Grid)[grid_index1];

					for(auto index_x = index_x_min; index_x <= index_x_max; index_x++)
						for(auto index_y = index_y_min; index_y <= index_y_max; index_y++)
							for(auto index_z = index_z_min; index_z <= index_z_max; index_z++) {

								PetscInt grid_index2 = FieldVar::GridID[i + index_x][j + index_y][k + index_z];

								if (grid_index2 >= 0) {
									PetscScalar tmp = .0;
									vector<PetscScalar> GridPos2 = (FieldVar::Grid)[grid_index2];

									for(auto atom = node->begin(); atom != node->end(); atom++) {

										PetscScalar aa_coords[] = {FieldVar::Coords_Local_x[*atom],
												           FieldVar::Coords_Local_y[*atom],
													   FieldVar::Coords_Local_z[*atom]};

										for(auto dim = 0; dim < FieldVar::Dim; dim++) {
											PetscScalar Jk = ((FieldVar::Grid)[grid_index2][dim] - aa_coords[dim]) * FieldVar::KernelFunction(aa_coords, FieldVar::Mass[*atom], GridPos2) * (2.0 / pow(FieldVar::Resol,2.0));
											tmp += Jk * ( (FieldVar::Grid)[grid_index1][dim] - aa_coords[dim]) / pow(FieldVar::Resol,2.0) * FieldVar::KernelFunction(aa_coords, FieldVar::Mass[*atom], GridPos) ;
										}
									}

									tmp *= 2.0;
									ierr = MatSetValues(FieldVar::Jacobian, 1, &grid_index1, 1, &grid_index2, &tmp, INSERT_VALUES);
								}
							}

						node++;
			}
		}

	MatAssemblyBegin(FieldVar::Jacobian, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(FieldVar::Jacobian, MAT_FINAL_ASSEMBLY);

        VecRestoreArray(Coords[0], &(FieldVar::Coords_Local_x));
        VecRestoreArray(Coords[1], &(FieldVar::Coords_Local_y));
        VecRestoreArray(Coords[2], &(FieldVar::Coords_Local_z));

	return ierr;
}

