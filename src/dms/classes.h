#ifndef DMS_CLASSES_H
#define DMS_CLASSES_H

#include "dms_headers.h"

#ifdef __cplusplus

class DmsBase;
//int DmsBase::getNumSS();
//int DmsBase::getNumSSglob();

typedef PetscErrorCode (*ptrMap)(CVec, CVec, DmsBase& Dbase, std::fstream&);
typedef PetscErrorCode (*ptrMapVelo)(CVec, DmsBase&, std::fstream&);
std::vector< std::vector<PetscInt> > dmsReadTop(char* fname);

/* The micro_state modifies only atomistic properties whereas the meso_state class modified
 * the CG properties. For this reason, only public inheritance is used to construct meso from
 * micro class, and all micro variables are private.
 */

class Micro_state {
public:

	Micro_state(const t_state* state, const t_mdatoms* tmdatoms, const gmx_mtop_t* top,
		    const t_inputrec*, PetscInt Dim, MPI_Comm, ptrMap, const int, const real, PetscInt, PetscInt, DmsBase*, char*, char*, rvec);
	Micro_state() {};

	explicit Micro_state(const Micro_state&);
	Micro_state & operator=(const Micro_state&);
	virtual ~Micro_state();

	PetscErrorCode Sync_DMS_fromMD(DmsBase*);
	PetscErrorCode Sync_MD_fromDMS(DmsBase*);

	PetscInt Get_DOF() const { return DOF; }
	PetscInt Get_DOF_local() const { return DOF_local; }
	PetscInt Get_Dim() const { return Dim; }

	CVec Get_Coords() const { return Coords; }
	//Cvec Get_cCoords() const { return cCoords; } Moved this to the son class Meso_State
	CVec Get_pCoords() const { return pCoords; }
	CVec Get_RefCoords() const { return Ref_Coords; }

	CVec Get_Velocities() const { return Velocities; }
	CVec Get_pVelocities() const { return pVelocities; }
	CVec Get_Forces() const { return Forces; }

	MPI_Comm Get_COMM() const { return COMM; }
	ptrMap mapping; // specific to the object: fine-graining, coarse-graining, etc.
	ptrMapVelo mappingVelo;

	const t_state* Get_MD_state() const { return MD_state; }
	const gmx_mtop_t* Get_MD_top() const { return MD_top; }
	Vec getMass() const { return Mass; }
	real getLength() const { return mdLength; }

	std::vector< std::vector<PetscInt> > getTopIndices() const { return topIndices; }
	Vec getEqLength() const { return eqLength; }

protected:
	MPI_Comm COMM; // the micro has its own comm because the atoms are parallelized uniquely

	CVec Coords,
	     pCoords,
	     cCoords, // used to save the constrained  
	     Velocities,
	     pVelocities,
	     Forces,
	     Ref_Coords;

	Vec Mass, eqLength;
	std::vector<PetscInt> atomIndices;
	PetscErrorCode ierr;

	PetscInt DOF,
		 Dim,
		 DOF_local, // this used to represent number of atoms per proc, but now is being used per subsystem
		 MPI_Size,
		 MPI_Rank,
		 istart, // to be set/used later??
		 iend;

private:
	// These data structures are stored on the master processor ONLY!!!
	const t_state* MD_state; // micro_state contains all the relevant info (natoms,
							 // positions, velocities, etc.) obtained from GROMACS
	const t_state* MD_state_local; // same as above but specific to each proc
	const t_mdatoms* mdatoms; // contains the local no of atoms on each proc ... among other things (??)
	const gmx_mtop_t* MD_top; // topology state
	const t_mdatoms* tmdatoms;
	real mdLength; // length of the MD micro phase in ps
	std::vector< std::vector<PetscInt> > topIndices;

	PetscErrorCode setupRefTop(char*, DmsBase*);
	PetscErrorCode setupSubsystem(char*);

	PetscInt *ssIndices;
	PetscInt numSS;
	PetscInt ssIndex;
	char* selFname;

	char* mode; // determines the mechanism by which a microstate is reconstructed.
	// Can be either default (Newtonian) or new (Lagrangian) mode
};


class Meso_state: public Micro_state {
	std::vector< CVec > cgTensor;
	PetscInt nHist;
public:
    
    Cvec Get_cCoords() const { return cCoords; }
    
	Meso_state(PetscInt NumCG, PetscInt DimCG, PetscInt, MPI_Comm, ptrMap, ptrMapVelo);

	Meso_state() {};

	explicit Meso_state(const Meso_state&);
	Meso_state &operator=(const Meso_state&);
	
	std::vector< CVec >& getCGTensor() { return cgTensor; }

	virtual ~Meso_state();
};

#endif
#endif
