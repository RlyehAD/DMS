#include "common.h"
#include "swm.h"
#include "petscdt.h"

/* Some terms/variables specific to the SWM:
 *
 * U: meso_micro_map of orthogonal polynomials
 * M: diagonal matrix of atomic masses
 * M * U: basis_weighted
 * ^t: transpose
 * r: coords (x,y,z)
 * (MU)^t * U : micro_meso_map
 *
 */

PetscErrorCode swm::constructBasis(const std::vector<Vec>& Coords, DmsBase& Dbase) {

	/* This function could be be called in parallel.
	 * The basis matrix is of size Natoms x nCG, and it is partitioned row-wise.
	 */
	PetscFunctionBegin;
	PetscErrorCode ierr;

	// construct reference config
	// very hackish!!!!
	ierr = Dbase.constructRef();
	CHKERRQ(ierr);

	std::vector<PetscScalar> box_size = Dbase.getBox(Coords);
	Mat *mesoMicroMap = Dbase.getMesoMicro();
	PetscInt istart, iend;

	ierr = MatGetOwnershipRange(*mesoMicroMap, &istart, &iend);
	CHKERRQ(ierr);

	std::vector<PetscInt> atomic_indices(iend - istart);
	std::vector< std::vector <PetscScalar> > polynomials(Dbase.Microscopic->Get_Dim());

	for(int i = 0; i < polynomials.size(); i++)
		polynomials[i].resize(iend - istart);

	for(int index = istart; index < iend; index++)
		atomic_indices[index - istart] = index;

	std::vector<Vec> NCoords(Dbase.Microscopic->Get_Dim());
	std::vector<PetscScalar> COM = Dbase.compCentOfMass(Coords);

	// The parallelization here must be checked thoroughly for efficiency

	for(int dim = 0; dim < Dbase.Microscopic->Get_Dim(); dim++) {

		// Compute U_{k1k2k3} s.t. k1 + k2 + k3 <= kmax
		ierr = VecCreateSeq(Dbase.getComm(), Dbase.Microscopic->Get_DOF_local(), NCoords.data() + dim);
		CHKERRQ(ierr);

		ierr = VecCopy(Coords[dim], NCoords[dim]);
		CHKERRQ(ierr);

		// Subtract center of gravity from coords
		PetscScalar COG;
		VecSum(NCoords[dim], &COG);
		COG /= Dbase.Microscopic->Get_DOF_local();

		ierr = VecShift(NCoords[dim], - COG);
		CHKERRQ(ierr);

		// Normalize coordinates to [-1,1]
		ierr = VecScale(NCoords[dim], 1.0/box_size[dim]);
		CHKERRQ(ierr);

		ierr = VecAssemblyBegin(NCoords[dim]);
		CHKERRQ(ierr);
		ierr = VecAssemblyEnd(NCoords[dim]);
		CHKERRQ(ierr);
	}

	PetscScalar *NCoords_ptr;
	std::vector< std::vector<PetscInt> > indices = computeIndices(Dbase.getOrder());

	for(int order = 0; order < Dbase.getNcg(); order++) {
		for(int dim = 0; dim < Dbase.Microscopic->Get_Dim(); dim++) {

			//VecView(NCoords[dim], PETSC_VIEWER_STDOUT_SELF);

			ierr = VecGetArray(NCoords[dim], &NCoords_ptr);
			CHKERRQ(ierr);

			ierr = PetscDTLegendreEval(iend - istart, NCoords_ptr, 1,
			       &indices[order][dim], polynomials[dim].data(), NULL, NULL);
			CHKERRQ(ierr);

			ierr = VecRestoreArray(NCoords[dim], &NCoords_ptr);
			CHKERRQ(ierr);
		}

		// Begin filling in the basis matrix
		for(int dim = 1; dim < Dbase.Microscopic->Get_Dim(); dim++)
				for(int i = 0; i < polynomials[0].size(); i++)
					polynomials[0][i] *= polynomials[dim][i];

		PetscScalar *poly_address = polynomials[0].data();

		ierr = MatSetValues(*mesoMicroMap, iend - istart, atomic_indices.data(), 1,
				    &order, poly_address, INSERT_VALUES);
		CHKERRQ(ierr);
	}

	ierr = MatAssemblyBegin(*mesoMicroMap, MAT_FINAL_ASSEMBLY);
	CHKERRQ(ierr);

	ierr = MatAssemblyEnd(*mesoMicroMap, MAT_FINAL_ASSEMBLY);
	CHKERRQ(ierr);

	/* This functions computed micro_meso_map = (MU)^t U.
	 * It is called whenever the reference configuration is updated.
	 */

	ierr = scaleBasis(Dbase);
	CHKERRQ(ierr);

	for(int dim = 0; dim < NCoords.size(); dim++) {
		ierr = VecDestroy(&NCoords[dim]);
		CHKERRQ(ierr);
	}

	PetscFunctionReturn(ierr);
}

PetscErrorCode swm::scaleBasis(DmsBase& Dbase) {
	PetscFunctionBegin;
	PetscErrorCode ierr;

	Mat *microMesoMap = Dbase.getMicroMeso(),
            *mesoMicroMap = Dbase.getMesoMicro(),
	    *kernel	  = Dbase.getKernel();

	ierr = MatCopy(*mesoMicroMap, *kernel, SAME_NONZERO_PATTERN);
	CHKERRQ(ierr);

	// Compute M x U
	ierr = MatDiagonalScale(*kernel, Dbase.Microscopic->getMass(), NULL);
	CHKERRQ(ierr);

	// Compute micro_meso_map => (MU)^T x U
	Mat MUT;
	ierr = MatTranspose(*kernel, MAT_INITIAL_MATRIX, &MUT);
	CHKERRQ(ierr);

	// Be careful with memory leaks here for micro_meso_map
	if(!(microMesoMap)) {
		ierr = MatDestroy(microMesoMap);
		CHKERRQ(ierr);
	}

	ierr = MatMatMult(MUT, *mesoMicroMap, MAT_INITIAL_MATRIX, PETSC_DEFAULT, microMesoMap);
	CHKERRQ(ierr);

	ierr = MatDestroy(&MUT);
	CHKERRQ(ierr);

	PetscFunctionReturn(ierr);
}

PetscErrorCode swm::coarseGrainVelo(CVec Coords, DmsBase& Dbase, std::fstream&)
{
        /* Solve (MU)^t * U \phi = (MU)^t * r
           i.e.  microMesoMap * Meso->Coords =  mesoMicroMap * Micro-Coords
         */

        PetscFunctionBegin;
        PetscErrorCode ierr;

        std::vector<PetscScalar> com = Dbase.compCentOfMass(Dbase.Microscopic->Get_RefCoords());

        Mat *microMesoMap = Dbase.getMicroMeso(),
            *mesoMicroMap = Dbase.getMesoMicro(),
            *kernel       = Dbase.getKernel();

        for(auto dim = 0; dim < Dbase.Mesoscopic->Get_Dim(); dim++) {
                Vec RHS;
                ierr = VecCreateSeq(Dbase.getComm(), Dbase.getNcg(), &RHS);
                CHKERRQ(ierr);

                KSP ksp;
                ierr = KSPCreate(Dbase.getComm(), &ksp);
                CHKERRQ(ierr);

                ierr = KSPSetOperators(ksp, *microMesoMap, *microMesoMap);
                CHKERRQ(ierr);

                ierr = KSPSetType(ksp, KSPCG);
                CHKERRQ(ierr);

                ierr = KSPSetFromOptions(ksp);
                CHKERRQ(ierr);

                ierr = KSPSetUp(ksp);
                CHKERRQ(ierr);

                ierr = MatMultTranspose(*kernel, Coords[dim], RHS);
                CHKERRQ(ierr);

                ierr = KSPSolve(ksp, RHS, Dbase.Mesoscopic->Get_Velocities()[dim]);
                CHKERRQ(ierr);

                // Free memory
                ierr = VecDestroy(&RHS);
                CHKERRQ(ierr);

                ierr = KSPDestroy(&ksp);
                CHKERRQ(ierr);
        }

        PetscFunctionReturn(ierr);
}

PetscErrorCode swm::coarseGrain(CVec Coords, CVec pCoords, DmsBase& Dbase, std::fstream&)
{
	/* Solve (MU)^t * U \phi = (MU)^t * r
	   i.e.  microMesoMap * Meso->Coords =  mesoMicroMap * Micro-Coords
	 */

	PetscFunctionBegin;
	PetscErrorCode ierr;

	std::vector<PetscScalar> com = Dbase.compCentOfMass(Dbase.Microscopic->Get_RefCoords());
 
	Mat *microMesoMap = Dbase.getMicroMeso(),
	    *mesoMicroMap = Dbase.getMesoMicro(),
	    *kernel	  = Dbase.getKernel();

	for(auto dim = 0; dim < Dbase.Mesoscopic->Get_Dim(); dim++) {
		Vec RHS;
		ierr = VecCreateSeq(Dbase.getComm(), Dbase.getNcg(), &RHS);
		CHKERRQ(ierr);

		KSP ksp;
		ierr = KSPCreate(Dbase.getComm(), &ksp);
		CHKERRQ(ierr);

		ierr = KSPSetOperators(ksp, *microMesoMap, *microMesoMap);
		CHKERRQ(ierr);

		ierr = KSPSetType(ksp, KSPCG);
		CHKERRQ(ierr);

		ierr = KSPSetFromOptions(ksp);
		CHKERRQ(ierr);

		ierr = KSPSetUp(ksp);
		CHKERRQ(ierr);

		ierr = VecShift(Coords[dim], -com[dim]);
		CHKERRQ(ierr);

		ierr = MatMultTranspose(*kernel, Coords[dim], RHS);
		CHKERRQ(ierr);

		ierr = KSPSolve(ksp, RHS, Dbase.Mesoscopic->Get_Coords()[dim]);
		CHKERRQ(ierr);

		//VecView(Dbase.Mesoscopic->Get_Coords()[dim], PETSC_VIEWER_STDOUT_SELF);

		ierr = VecShift(Coords[dim], com[dim]);
                CHKERRQ(ierr);

		// Free memory
		ierr = VecDestroy(&RHS);
		CHKERRQ(ierr);

		ierr = KSPDestroy(&ksp);
		CHKERRQ(ierr);
	}

	PetscFunctionReturn(ierr);
}
