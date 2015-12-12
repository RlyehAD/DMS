											/********************************************************/
											/****** Continuum Field Variables: Coarse-graining *****/
											/******************************************************/

/* The idea here is to coarse grain an atomistic state (the fine scale state) by defining a kernel function K such that:
 *
 * 										FV = \int K * FV / density
 *
 * The choice of the field variable (FV) used here is the density i.e. FV becomes in discretized form:
 *
 * 											FV = \sum_i K_i
 *
 * Based on this the FV momenta and forces are easily computed.
 *
 * TODO: This file ought to be optimized by combining the force,vel, pos for FVs evaluations into a single function that avoids three time loops over
 * the AtomicList.
 *
 * TODO: ??? worth parallelizing or stick to serial code ???
 */


#include "Interface.h"

PetscScalar FieldVar::KernelFunction(const PetscScalar* const Coords, PetscInt mass, const vector<PetscScalar> &GridPos) {
	PetscScalar dotProd = .0;

	for(int d = 0; d < GridPos.size(); d++)
		dotProd += pow(GridPos[d] - Coords[d], 2.0);

	dotProd = dotProd / pow(FieldVar::Resol, 2.0);

	return mass * exp(-dotProd);
}

vector<PetscScalar> FieldVar::ComputeFV(const PetscScalar* const Coords) {

	/***********************************************************************/
	/****** Computes the field variables of choice: MPI not YET DONE! *****/
	/*********************************************************************/

	PetscFunctionBegin;
	auto it = FieldVar::FVList.begin(); // FVList contains atomic indices for each FV cell / node

	for(auto i = 0; i < FieldVar::AdjNumNodes; i++) {
		FieldVar::FieldVars[i] = 0.f; // clear previously stored values

		vector<PetscScalar> GridPos = FieldVar::Grid[i];

		for(auto aa_it = it->begin(); aa_it != it->end(); aa_it++)
			FieldVar::FieldVars[i] += FieldVar::KernelFunction(Coords+(FieldVar::Dim)*(*aa_it), FieldVar::Mass[*aa_it], GridPos);

		it++;
	}
	return FieldVar::FieldVars;
}

vector<PetscScalar> FieldVar::ComputeFV_Vel(const PetscScalar* const Coords, const PetscScalar* const Vel) {

	/*************************************************************************/
	/****** Computes the velocity of field variables: MPI not YET DONE! *****/
	/***********************************************************************/

	PetscFunctionBegin;
	auto it = FieldVar::FVList.begin();

	for(auto i = 0; i < FieldVar::AdjNumNodes; i++) {

		FieldVar::FieldVars_Vel[i] = 0.f;
		vector<PetscScalar> GridPos = FieldVar::Grid[i]; // must be deleted when re-equilibriating to avoid memory leaks

		for(auto aa_it = it->begin(); aa_it != it->end(); aa_it++) {

			const PetscScalar *const AtomicCoords = Coords+(FieldVar::Dim)*(*aa_it),
					  *const AtomicVel = Vel+(FieldVar::Dim)*(*aa_it),
					  Kernel = FieldVar::KernelFunction(AtomicCoords, FieldVar::Mass[*aa_it], GridPos);

			for(auto d = 0; d < FieldVar::Dim; d++)
				FieldVar::FieldVars_Vel[i] += 2.0 * (GridPos[d] - AtomicCoords[d]) / pow(FieldVar::Resol,2.0) * AtomicVel[d] * Kernel;
		}

		it++;
	}

	return FieldVar::FieldVars_Vel;
}

vector<PetscScalar> FieldVar::ComputeFV_For(const PetscScalar* const Coords, const PetscScalar* const Vel, const PetscScalar* const For) {

	/**********************************************************************/
	/****** Computes the force of field variables: MPI not YET DONE! *****/
	/********************************************************************/

	PetscFunctionBegin;
	auto it = FieldVar::FVList.begin();

	for(auto i = 0; i < FieldVar::AdjNumNodes; i++) {

		FieldVar::FieldVars_For[i] = 0.f;
		vector<PetscScalar> GridPos = FieldVar::Grid[i]; // must be deleted when re-equilibriating to avoid memory leaks

		for(auto aa_it = it->begin(); aa_it != it->end(); aa_it++) {

			const PetscScalar *const AtomicCoords = Coords+(FieldVar::Dim)*(*aa_it),
					  *const AtomicVel = Vel+(FieldVar::Dim)*(*aa_it),
					  *const AtomicFor = For+(FieldVar::Dim)*(*aa_it),
					  Kernel = FieldVar::KernelFunction(AtomicCoords, FieldVar::Mass[*aa_it], GridPos);

			for(auto d = 0; d < FieldVar::Dim; d++)
				FieldVar::FieldVars_For[i] += 2.0 * ( (GridPos[d] - AtomicCoords[d]) / pow(FieldVar::Resol,2.0) * AtomicFor[d]
					- AtomicVel[d] / pow(FieldVar::Resol,2.0) * AtomicVel[d]
					+ 2.0 * AtomicVel[d] * pow((GridPos[d] - AtomicCoords[d]) / pow(FieldVar::Resol,2.0), 2.0) * AtomicVel[d] ) * Kernel;
		}

		it++;
	}

	return FieldVar::FieldVars_For;
}
