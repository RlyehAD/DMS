#include "../dms_headers.h"

#undef __FUNCT__
#define __FUNCT__ "ComputeMetric"
PetscErrorCode ComputeMetric(std::vector<Vec> Coords, std::vector< std::vector<PetscInt> > Indices, Vec& EqLength, Vec& Constraints) {
  /* ** Description **
     This function computes an array (distriuted on all procs) of equations satisfying
     the difference between the metric squared of two points in space and their
     equilibrium metric squared. This metric in this case corresponds to atomic
     distance (that is, for bonds/angles).

     |r_i - r_j|^2.0 - (l_ij)^2.0 with l_ij being the equilibrium metric to be satisfied.

     ** Credit **
     Andrew Abi Mansour
     Department of Chemistry
     Indiana University, Bloomington

     ** Last update **
     May 11, 2013 - 23:15
  */
	PetscFunctionBegin;

	std::cout << __FUNCT__ << std::endl;
	PetscErrorCode ierr;

	// Extract appropriate dimension
	PetscInt NumAtoms, Dims;
	Dims = Coords.size();
	VecGetSize(Coords[0], &NumAtoms);

	// Get local partition of the Constraints vector
	int istart, iend;
	VecGetOwnershipRange(Constraints,&istart,&iend);

	PetscInt    *index         = new PetscInt[iend - istart];
	PetscScalar *atomCoordsOne = new PetscScalar[iend - istart],
	            *atomCoordsTwo = new PetscScalar[iend - istart];

	Vec tmp;

	ierr = VecZeroEntries(Constraints);
	CHKERRQ(ierr);

	for(auto dim = 0; dim < Dims; dim++) {

		for(auto i = 0; i < iend - istart; i++)
                	index[i] = Indices[i][0];

		ierr = VecGetValues(Coords[dim], iend - istart, index, atomCoordsOne); 
		CHKERRQ(ierr);

		for(auto i = 0; i < iend - istart; i++)  
                        index[i] = Indices[i][1];

		ierr = VecGetValues(Coords[dim], iend - istart, index, atomCoordsTwo); CHKERRQ(ierr);

		for(auto i = 0; i < iend - istart; i++)
			atomCoordsTwo[i] = pow(atomCoordsOne[i] - atomCoordsTwo[i], 2.0);

		ierr = VecCreateSeqWithArray(PETSC_COMM_SELF, 1, iend - istart, atomCoordsTwo, &tmp);
		ierr = VecAXPY(Constraints, 1.0, tmp);
		ierr = VecDestroy(&tmp);
	}

	ierr = VecAXPBY(Constraints, -1.0, 1.0, EqLength); // Cons = dx^y + dy^2 + dz^2 - Leq^2

	VecAssemblyBegin(Constraints);
	VecAssemblyEnd(Constraints);

	delete[] index;
	delete[] atomCoordsOne;
	delete[] atomCoordsTwo;

	PetscFunctionReturn(ierr);
}

#undef __FUNCT__
#define __FUNCT__ "ComputeAtomicJacobi"
PetscErrorCode ComputeAtomicJacobi(std::vector<Vec>  Coords, std::vector< std::vector<PetscInt> > Indices, std::vector<Mat*> AtomicJacobi) {
  /* ** Description **
        This function constructs the Jacobian of the metric equations.

	** Size **
	The matrix should be of size NumCons x 3*NumAtoms

	** Last updated **

	** Credit **
	Andrew Abi Mansour
    Department of Chemistry
    Indiana University, Bloomington
  */
	PetscFunctionBegin;

	std::cout << __FUNCT__ << std::endl;
	
	PetscErrorCode ierr;
	int istart, iend;
	int Dims, NumAtoms;

	ierr = MatGetOwnershipRange(*AtomicJacobi[0], &istart, &iend); CHKERRQ(ierr);

	Dims = Coords.size();
	VecGetSize(Coords[0], &NumAtoms);

	PetscInt    *index         = new PetscInt[iend - istart];
    	PetscScalar *atomCoordsOne = new PetscScalar[iend - istart],
            	    *atomCoordsTwo = new PetscScalar[iend - istart];

	PetscScalar vals[2];

	for(auto dim = 0; dim < Dims; dim++) {

		for(auto i = 0; i < iend - istart; i++)
			index[i] = Indices[i][0];

		ierr = VecGetValues(Coords[dim], iend - istart, index, atomCoordsOne); 
		CHKERRQ(ierr);

		for(auto i = 0; i < iend - istart; i++)
                        index[i] = Indices[i][1];

		ierr = VecGetValues(Coords[dim], iend - istart, index, atomCoordsTwo); 
		CHKERRQ(ierr);

		for(auto i = istart; i < iend; i++) {
			vals[0] = 2.0 * (atomCoordsOne[i] - atomCoordsTwo[i]);
			vals[1] = - 2.0 * (atomCoordsOne[i] - atomCoordsTwo[i]);

			ierr = MatSetValues(*AtomicJacobi[dim], 1, &i, 2, (Indices.data() + i)->data(), vals, INSERT_VALUES); 
			CHKERRQ(ierr);
		}

		ierr = MatAssemblyBegin(*AtomicJacobi[dim], MAT_FINAL_ASSEMBLY); 
		CHKERRQ(ierr);

  		ierr = MatAssemblyEnd(*AtomicJacobi[dim], MAT_FINAL_ASSEMBLY); 
		CHKERRQ(ierr);
	}

  	delete[] index;
  	delete[] atomCoordsOne;
  	delete[] atomCoordsTwo;

	PetscFunctionReturn(ierr);
}

#undef __FUNCT__
#define __FUNCT__ "ComputeAdjacencyMatrix"
PetscErrorCode ComputeAdjacencyMatrix(std::vector< std::vector<PetscInt> > Indices, Mat& Adj, const PetscInt Dims) {
  /* 
     ** Description ** 
     This function computes an adjacency-like matrix which has the entries 1 and -1
     for each row: +1 for the first metric index, -1 for the second, i.e.
     
     \sum_d=1^3 Diag(P x r_d) x (P x r_d) = Metric^2.0

     ** Size **
     The total size of Adj should be NumCons x NumAtoms.

     ** Technical details **
     In the spirit of communication-free algorithms, Adj is assembled locally on each
     processor. For efficient storage and floating-point operations, the matrix is
     stored in sparse format. 

     ** Last updated **
     May 10, 2013 - 11:45

     ** Credit **
     Andrew Abi Mansour
     Department of Chemistry
     Indiana University, Bloomington
  */
	PetscFunctionBegin;

	std::cout << __FUNCT__ << std::endl;

	PetscErrorCode ierr;
	PetscInt NumCons, NumAtoms, istart, iend;

	// Get the sizes of Indices (NumCons x 2)
	ierr = MatGetSize(Adj, &NumCons, &NumAtoms);  CHKERRQ(ierr);
	ierr = MatGetOwnershipRange(Adj,&istart,&iend);  CHKERRQ(ierr);

	PetscScalar vals[] = {1.0, -1.0};

	for(auto i = istart; i < iend; i++)
		ierr = MatSetValues(Adj, 1, &i, 2, (Indices.data() + i)->data(), vals, INSERT_VALUES); CHKERRQ(ierr);

  	MatAssemblyBegin(Adj, MAT_FINAL_ASSEMBLY);
  	MatAssemblyEnd(Adj, MAT_FINAL_ASSEMBLY);

  	PetscFunctionReturn(ierr);
}

#undef __FUNCT__
#define __FUNCT__ "AssembleLagJacobi"
PetscErrorCode AssembleLagJacobi(std::vector<Vec> Coords, Mat& Adjacency, Mat& AdjacencyTrans, std::vector<Mat*> AtomicJacobi, Mat& LagJacobiGlobal) {
       PetscFunctionBegin;

       std::cout << __FUNCT__ << std::endl;

       PetscErrorCode ierr;
       PetscInt Dims, NumAtoms, NumCons, istart, iend;
       MatGetSize(AdjacencyTrans, &NumAtoms, &NumCons);

       VecGetSize(Coords[0], &NumAtoms);
       Dims = Coords.size();

       MatGetOwnershipRange(Adjacency, &istart, &iend);

       for(auto dim = 0; dim < Dims; dim++) {

	       // Compute Adj * Coords product
	       Vec Ar;
	       ierr = VecCreateSeq(PETSC_COMM_SELF, NumCons, &Ar); CHKERRQ(ierr);
	       ierr = MatMult(Adjacency, Coords[dim], Ar); CHKERRQ(ierr);

	       ierr = VecScale(Ar, 2.0);

	       // Compute Jc * Adj.T product
	       Mat AdjAtom;
	       ierr = MatCreate(PETSC_COMM_SELF, &AdjAtom); CHKERRQ(ierr);
	       ierr = MatSetSizes(AdjAtom, PETSC_DECIDE, PETSC_DECIDE, NumCons, NumCons); CHKERRQ(ierr);
	       ierr = MatSetType(AdjAtom, MATMPIAIJ); CHKERRQ(ierr);

	       ierr = MatMatMult(*AtomicJacobi[dim], AdjacencyTrans, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &AdjAtom); CHKERRQ(ierr);

	       MatDiagonalScale(AdjAtom, NULL, Ar);

	       if(dim == 0) {

		    MatDuplicate(AdjAtom, MAT_COPY_VALUES ,&LagJacobiGlobal);
	       	    CHKERRQ(ierr);
	       }

	       else {

	            ierr = MatAXPY(LagJacobiGlobal, 1.0, AdjAtom, SAME_NONZERO_PATTERN);
	       	    CHKERRQ(ierr);
	       }

	       // Destroy allocated vectors & matrices
	       VecDestroy(&Ar);
	       MatDestroy(&AdjAtom);
	   }

       return ierr;
}

#undef __FUNCT__
#define __FUNCT__ "WriteMatrixToFile"
void WriteMatrixToFile(Mat& Matrix, char* ofname) {
        PetscFunctionBegin;
	PetscViewer viewer;
	PetscViewerASCIIOpen(PETSC_COMM_SELF, ofname, &viewer);
	PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
	MatView(Matrix, viewer);
	PetscViewerDestroy(&viewer);
}

#undef __FUNCT__
#define __FUNCT__ "consOptim"
PetscErrorCode consOptim(std::vector<Vec> Coords, std::vector<std::vector< PetscInt> >& Indices, Vec EqLength, float scale) {
     /* 
        ** Description **
	This is the main computational engine that solves the equations

	\nabla_r f = r - r_u - Jc.T * L = 0,
	\nabla_L f = Diag(Adj * r) Adj * r - l_eq^2.0 = 0.
	
	** Credit **
	Andrew Abi Mansour
	Department of Chemistry
	Indiana University, Bloomington
     */

	PetscFunctionBegin;

	PetscErrorCode ierr;
	PetscReal tol1 = 0.01, tol2 = 0.1, error1 = tol1, error2 = tol2;
	PetscInt NumCons, NumAtoms, Dims, ConsDims = 2, rank, size, istart, iend;

	MPI_Comm_size(PETSC_COMM_WORLD, &size);
  	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	Dims = Coords.size();
	VecGetSize(Coords[0], &NumAtoms);
	VecGetSize(EqLength, &NumCons);

  	Vec Constraints, mu, dr;
 	VecCreateMPI(PETSC_COMM_SELF, PETSC_DECIDE, NumCons, &Constraints);
 	VecCreateMPI(PETSC_COMM_SELF, PETSC_DECIDE, NumCons, &mu);
 	VecCreateMPI(PETSC_COMM_SELF, PETSC_DECIDE, NumAtoms, &dr);

	std::cout << __FUNCT__ << std::endl;

	// Create the Seidel adjacency matrix locally (on each processor)
	Mat Adjacency;
	MatCreate(PETSC_COMM_SELF, &Adjacency);
	MatSetSizes(Adjacency, PETSC_DECIDE, PETSC_DECIDE, NumCons, NumAtoms);
	MatSetType(Adjacency, MATMPIAIJ);
	MatMPIAIJSetPreallocation(Adjacency, ConsDims, PETSC_NULL, ConsDims, PETSC_NULL);
	ComputeAdjacencyMatrix(Indices, Adjacency, Dims);

 	// Create AtomicJacobi and distribute it across all processors

	Mat AtomicJacobiX, AtomicJacobiY, AtomicJacobiZ;
	MatCreate(PETSC_COMM_SELF, &AtomicJacobiX);
	MatSetSizes(AtomicJacobiX, PETSC_DECIDE, PETSC_DECIDE, NumCons, NumAtoms);
	MatSetType(AtomicJacobiX, MATMPIAIJ);
	MatMPIAIJSetPreallocation(AtomicJacobiX, ConsDims, PETSC_NULL, ConsDims, PETSC_NULL);

        MatCreate(PETSC_COMM_SELF, &AtomicJacobiY);
        MatSetSizes(AtomicJacobiY, PETSC_DECIDE, PETSC_DECIDE, NumCons, NumAtoms);
        MatSetType(AtomicJacobiY, MATMPIAIJ);
        MatMPIAIJSetPreallocation(AtomicJacobiY, ConsDims, PETSC_NULL, ConsDims, PETSC_NULL);

        MatCreate(PETSC_COMM_SELF, &AtomicJacobiZ);
        MatSetSizes(AtomicJacobiZ, PETSC_DECIDE, PETSC_DECIDE, NumCons, NumAtoms);
        MatSetType(AtomicJacobiZ, MATMPIAIJ);
        MatMPIAIJSetPreallocation(AtomicJacobiZ, ConsDims, PETSC_NULL, ConsDims, PETSC_NULL);

	std::vector<Mat*> AtomicJacobi(Dims);
	AtomicJacobi[0] = &AtomicJacobiX;
	AtomicJacobi[1]	= &AtomicJacobiY;
	AtomicJacobi[2]	= &AtomicJacobiZ;

 	// Create linear solver
 	KSP ksp;
 	KSPCreate(PETSC_COMM_SELF, &ksp);
	KSPSetInitialGuessNonzero(ksp,PETSC_FALSE);
	KSPSetType(ksp,KSPIBCGS);

	//KSPSetTolerances(ksp,0.01,0.001,0.1,1000);
 	VecGetOwnershipRange(dr,&istart,&iend);

	Mat AdjacencyTrans;
	MatCreate(PETSC_COMM_SELF, &AdjacencyTrans);
	MatSetSizes(AdjacencyTrans, PETSC_DECIDE, PETSC_DECIDE, NumAtoms, NumCons);
	MatSetType(AdjacencyTrans, MATMPIAIJ);
	MatTranspose(Adjacency, MAT_INITIAL_MATRIX, &AdjacencyTrans);

 	while(error2 >= tol2) {
 		
 		ierr = ComputeMetric(Coords, Indices, EqLength, Constraints); CHKERRQ(ierr);

 		ierr = ComputeAtomicJacobi(Coords, Indices, AtomicJacobi); CHKERRQ(ierr);
	       
 		Mat LagJacobi;
		ierr = AssembleLagJacobi(Coords, Adjacency, AdjacencyTrans, AtomicJacobi, LagJacobi); CHKERRQ(ierr);
		ierr = MatShift(LagJacobi, scale); //NumCons / 100.0);
		CHKERRQ(ierr);
		
		// The nnz pattern does not change because the macromolecule
		// topology does not change.
		
 		ierr = KSPSetOperators(ksp,LagJacobi,LagJacobi); CHKERRQ(ierr);
 		ierr = KSPSolve(ksp,Constraints,mu); CHKERRQ(ierr);
		
 		// Update r = ru + dr where dr = - Jt * mu

 		for(auto dim = 0; dim < Dims; dim++) {

 			ierr = MatMultTranspose(*AtomicJacobi[dim], mu, dr); 
			CHKERRQ(ierr);

			ierr = VecAXPBY(Coords[dim], -1.0, 1.0, dr);
			CHKERRQ(ierr);

 			ierr = VecAssemblyBegin(Coords[dim]); CHKERRQ(ierr);
  			ierr = VecAssemblyEnd(Coords[dim]); CHKERRQ(ierr);
  		}

 		VecNorm(dr, NORM_INFINITY, &error1);
 		VecNorm(Constraints, NORM_INFINITY, &error2); 
		MatDestroy(&LagJacobi);
 		std::cout << "atom disp = " << error1 <<  " " << tol1 << " " << error2 << " " << tol2 << std::endl;
 	}

 	//MatView(LagJacobi, PETSC_VIEWER_STDOUT_WORLD);
  	//VecView(Constraints,PETSC_VIEWER_STDOUT_SELF);
	
	for(auto dim = 0; dim < Dims; dim++)
  		MatDestroy(AtomicJacobi[dim]);

	MatDestroy(&AdjacencyTrans);
	MatDestroy(&Adjacency);

 	VecDestroy(&Constraints);
	VecDestroy(&dr);
	VecDestroy(&mu);

	PetscFunctionReturn(ierr);
}
