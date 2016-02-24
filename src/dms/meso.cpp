#include "classes.h"

Meso_state::Meso_state(PetscInt NumCG, PetscInt DimCG, PetscInt nHistory, MPI_Comm Comm, ptrMap ptrFunc, ptrMapVelo ptrFuncVelo) {

	PetscFunctionBegin;

	DOF = NumCG;
	Dim = DimCG;
	COMM = Comm;
	mapping = ptrFunc;
	mappingVelo = ptrFuncVelo;
	nHist = nHistory;

	/*
	cgTensor.resize(nHist);

	for(int fr = 0; fr < nHist; fr++) {
		cgTensor[fr].resize(Dim);
		
		for(int dim = 0; dim < Dim; dim++)
			VecCreateMPI(COMM, PETSC_DECIDE, DOF, &cgTensor[fr][dim]);
 	}
	*/

	Coords.resize(Dim);
	Ref_Coords.resize(Dim);
	pCoords.resize(Dim);
	Velocities.resize(Dim);
	pVelocities.resize(Dim);
	Forces.resize(Dim);

	for(int dim = 0; dim < Dim; dim++)

		if(COMM == PETSC_COMM_WORLD) {
			VecCreateMPI(COMM, PETSC_DECIDE, NumCG, Coords.data() + dim);
			VecZeroEntries(Get_Coords()[dim]);
			VecAssemblyBegin(Get_Coords()[dim]);
			VecAssemblyEnd(Get_Coords()[dim]);

			VecCreateMPI(COMM, PETSC_DECIDE, NumCG, pCoords.data() + dim);
			VecZeroEntries(Get_pCoords()[dim]);
			VecAssemblyBegin(Get_pCoords()[dim]);
			VecAssemblyEnd(Get_pCoords()[dim]);

			VecCreateMPI(COMM, PETSC_DECIDE, NumCG, Ref_Coords.data() + dim);
                        VecZeroEntries(Get_RefCoords()[dim]);
                        VecAssemblyBegin(Get_RefCoords()[dim]);
                        VecAssemblyEnd(Get_RefCoords()[dim]);

			VecCreateMPI(COMM, PETSC_DECIDE, NumCG, Velocities.data() + dim);
			VecZeroEntries(Get_Velocities()[dim]);
			VecAssemblyBegin(Get_Velocities()[dim]);
			VecAssemblyEnd(Get_Velocities()[dim]);

			VecCreateMPI(COMM, PETSC_DECIDE, NumCG, pVelocities.data() + dim);
                        VecZeroEntries(Get_pVelocities()[dim]);
                        VecAssemblyBegin(Get_pVelocities()[dim]);
                        VecAssemblyEnd(Get_pVelocities()[dim]);

			VecCreateMPI(COMM, PETSC_DECIDE, NumCG, Forces.data() + dim);
			VecZeroEntries(Get_Forces()[dim]);
			VecAssemblyBegin(Get_Forces()[dim]);
			VecAssemblyEnd(Get_Forces()[dim]);
		}
		else {
			VecCreateSeq(COMM, NumCG, Coords.data() + dim);
			VecZeroEntries(Get_Coords()[dim]);

			VecCreateSeq(COMM, NumCG, pCoords.data() + dim);
			VecZeroEntries(Get_pCoords()[dim]);

			VecCreateSeq(COMM, NumCG, Ref_Coords.data() + dim);
                        VecZeroEntries(Get_RefCoords()[dim]);

			VecCreateSeq(COMM, NumCG, Velocities.data() + dim);
			VecZeroEntries(Get_Velocities()[dim]);

			VecCreateSeq(COMM, NumCG, pVelocities.data() + dim);
                        VecZeroEntries(Get_pVelocities()[dim]);

			VecCreateSeq(COMM, NumCG, Forces.data() + dim);
			VecZeroEntries(Get_Forces()[dim]);
		}
}

Meso_state::Meso_state(const Meso_state& meso_class) {

	PetscFunctionBegin;

	DOF = meso_class.DOF;
	Dim = meso_class.Dim;
	COMM = meso_class.COMM;
	mapping = meso_class.mapping;

	Coords.resize(Dim);
	pCoords.resize(Dim);
	Velocities.resize(Dim);
	pVelocities.resize(Dim);
	Forces.resize(Dim);

	for(int dim = 0; dim < Dim; dim++)
		if(COMM == PETSC_COMM_WORLD) {
			VecCreateMPI(COMM, PETSC_DECIDE, DOF, Coords.data() + dim);
			VecCopy(meso_class.Get_Coords()[dim], Get_Coords()[dim]);
			VecAssemblyBegin(Get_Coords()[dim]);
			VecAssemblyEnd(Get_Coords()[dim]);

			VecCreateMPI(COMM, PETSC_DECIDE, DOF, pCoords.data() + dim);
			VecCopy(meso_class.Get_pCoords()[dim], Get_pCoords()[dim]);
			VecAssemblyBegin(Get_pCoords()[dim]);
			VecAssemblyEnd(Get_pCoords()[dim]);

			VecCreateMPI(COMM, PETSC_DECIDE, DOF, Velocities.data() + dim);
			VecCopy(meso_class.Get_Velocities()[dim], Get_Velocities()[dim]);
			VecAssemblyBegin(Get_Velocities()[dim]);
			VecAssemblyEnd(Get_Velocities()[dim]);

			VecCreateMPI(COMM, PETSC_DECIDE, DOF, pVelocities.data() + dim);
                        VecCopy(meso_class.Get_pVelocities()[dim], Get_pVelocities()[dim]);
                        VecAssemblyBegin(Get_pVelocities()[dim]);
                        VecAssemblyEnd(Get_pVelocities()[dim]);

			VecCreateMPI(COMM, PETSC_DECIDE, DOF, Forces.data() + dim);
			VecCopy(meso_class.Get_Forces()[dim], Get_Forces()[dim]);
			VecAssemblyBegin(Get_Forces()[dim]);
			VecAssemblyEnd(Get_Forces()[dim]);
		}
		else {
			VecCreateSeq(COMM, DOF, Coords.data() + dim);
			VecCopy(meso_class.Get_Coords()[dim], Get_Coords()[dim]);

			VecCreateSeq(COMM, DOF, pCoords.data() + dim);
			VecCopy(meso_class.Get_pCoords()[dim], Get_pCoords()[dim]);

			VecCreateSeq(COMM, DOF, Velocities.data() + dim);
			VecCopy(meso_class.Get_Velocities()[dim], Get_Velocities()[dim]);

			VecCreateSeq(COMM, DOF, pVelocities.data() + dim);
                        VecCopy(meso_class.Get_pVelocities()[dim], Get_pVelocities()[dim]);

			VecCreateSeq(COMM, DOF, Forces.data() + dim);
			VecCopy(meso_class.Get_Forces()[dim], Get_Forces()[dim]);
		}
}

Meso_state& Meso_state::operator=(const Meso_state& meso_class) {

	PetscFunctionBegin;

	COMM = meso_class.COMM;
	DOF = meso_class.DOF;
	Dim = meso_class.Dim;
	mapping = meso_class.mapping;

	Coords.resize(Dim);
	Velocities.resize(Dim);
	Forces.resize(Dim);

	// TODO: See if vectors have already been allocated memory
	//	 	 support for sequential computations
	//		 support for pCoords

	for(int dim = 0; dim < Dim; dim++) {

			VecCreateMPI(COMM, PETSC_DECIDE, DOF, Coords.data() + dim);
			VecCopy(meso_class.Get_Coords()[dim], Get_Coords()[dim]);
			VecAssemblyBegin(Get_Coords()[dim]);
			VecAssemblyEnd(Get_Coords()[dim]);

			VecCreateMPI(COMM, PETSC_DECIDE, DOF, Velocities.data() + dim);
			VecCopy(meso_class.Get_Velocities()[dim], Get_Velocities()[dim]);
			VecAssemblyBegin(Get_Velocities()[dim]);
			VecAssemblyEnd(Get_Velocities()[dim]);

			VecCreateMPI(COMM, PETSC_DECIDE, DOF, Forces.data() + dim);
			VecCopy(meso_class.Get_Forces()[dim], Get_Forces()[dim]);
			VecAssemblyBegin(Get_Forces()[dim]);
			VecAssemblyEnd(Get_Forces()[dim]);

		}

	// Must equate vectors here
	/*
	for(auto dim = 0; dim < DimCG; dim++) {
		Meso_state::Coords[dim] = meso_class.Coords[dim];
		Meso_state::Velocities[dim] = meso_class.Velocities[dim];
		Meso_state::Forces[dim] = meso_class.Forces[dim];
	}
	*/

	return *this;
}

Meso_state::~Meso_state() {
	PetscFunctionBegin;
}
