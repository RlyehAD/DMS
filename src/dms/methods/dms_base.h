 /*
 * This file is part of the DMS molecular simulation package.
 *
 * Copyright (c) 2014, by the DMS development team:,
 * Andrew Abi Mansour
 *
 * DMs is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * DMS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with DMS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund DMS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.dms.org

 @Author: Andrew Abi Mansour
 @Created: March 27, 2014

 *This file should include ONLY the functions to be interfaced between DMS (CG phase) and GROMACS (MD phase)
 *It should contain all the necessary functions related to the CG method of choice.
 *In order to interface classes with C, a void DMS pointer is created to handle
 *all the possible methods to be called from md.c

 *The DMS class provides the base class from which a specific CG method inherits
 *its properties.
*/

#ifndef CG_METHODS_H
#define CG_METHODS_H

typedef void* dmsBasePtr;

#ifdef __cplusplus

#include "../dms_headers.h"
#include "space_warping/swm.h"
#include "../classes.h"
#include "../integrators/factorization/euler/euler_integrator.h"
#include "../integrators/factorization/pade/pade_integrator.h"

Mat DMS_Mat_Create(const PetscInt nrows, const PetscInt ncols);
PetscErrorCode consOptim(std::vector<Vec> Coords, std::vector< std::vector<PetscInt> >&, Vec EqLength, float);

class DmsBase {
public:
	DmsBase(const t_state* state, const t_mdatoms* mdatoms,
                const gmx_mtop_t* top, const t_inputrec* ir, const gmx_int64_t aDim, const gmx_int64_t cDim, const int max_order,
                const gmx_int64_t freq, const real dt, const gmx_int64_t t0, MPI_Comm comm, const int mSteps, const double optimScale, const PetscInt, 
		const PetscInt, const PetscInt, std::string = "SpaceWarping", char* readref = NULL, char* topFname = NULL, char* subFname = NULL, rvec forces[] = NULL);

	DmsBase(const t_state* state, const t_mdatoms* tmdatoms,
			const gmx_mtop_t* top, const t_inputrec* ir, const gmx_int64_t aDim, const gmx_int64_t cDim, const gmx_int64_t nCGx,
			const gmx_int64_t nCGy, const gmx_int64_t nCGz, const gmx_int64_t freq, const gmx_int64_t assFreq, const real dt, MPI_Comm communic, 
			const int mSteps, const real resol, const real scaling, const gmx_int64_t thresh, const int, const PetscInt,
			const PetscInt ssIndex, std::string cgType = "FieldVariables", bool deb = false, char* topFname = NULL, char* subFname = NULL);

	Micro_state *Microscopic;
    	Meso_state  *Mesoscopic;

	virtual PetscInt cgStep(gmx_int64_t); // do we really need this to be virtual? Perhaps it's best to 
	// prevent overriding this in derived classes

	DmsIntegrator *Integrator;
	virtual ~DmsBase();

	bool IsDebug() { return debug; }
	const t_inputrec* tprRecord;
	void registerMethods();
	static std::string getTime();

	PetscErrorCode constructRef();
	PetscErrorCode constructVelocities(); // useless?
    	PetscErrorCode constructCoords();
	PetscErrorCode constructForces();
    
	// Getters
	gmx_int64_t getNatoms() const { return nAtoms; }
	gmx_int64_t getLocalNatoms() const { return nLocalAtoms; }
	gmx_int64_t getNcg() const { return nCG; }
	gmx_int64_t getNcgAdj() const { return nCGAdj; }
	gmx_int64_t getNcgX() const { return nCGx; }
	gmx_int64_t getNcgY() const { return nCGy; }
	gmx_int64_t getNcgZ() const { return nCGz; }
	MPI_Comm getComm() const { return comm; }
	PetscInt getRank() const { return mpiRank; }
	gmx_int64_t getTimeStep() const { return timeStep; }
	gmx_int64_t getFreqUpdate() const { return freqUpdate; }	
    	PetscInt getOrder() const { return order; }
	gmx_int64_t getThresh() const { return threshold; }
	gmx_int64_t getNeighX() const { return neighX; }
	gmx_int64_t getNeighY() const { return neighY; }
	gmx_int64_t getNeighZ() const { return neighZ; }
	gmx_int64_t getAssemFreq() const { return assembleFreq; }
	std::vector<PetscScalar> compCentOfMass(const CVec& Coords);
	std::vector<PetscScalar> getBox(const CVec& Coords) const;

	int getNumSS();
	int getNumSSglob(); 

	std::vector<Mat>& getMesoMicroVec() { return mesoMicroVec; }
	std::vector<Mat>& getMicroMesoVec() { return microMesoVec; }

        Mat* getMesoMicro() { return &mesoMicroMap; }
        Mat* getMicroMeso() { return &microMesoMap; }

	Mat* getKernel() { return &kernel; }
	Mat* getKernelVec(int n) { return &kernelVec[n]; }

	std::vector< std::vector<PetscScalar> >& getGrid() { return grid; }
	std::vector< std::list<PetscInt> >& getFVlist() { return fvList; }
	std::vector< std::vector< std::vector<PetscInt> > >& getGridID() { return gridID; }

	// Setters
	void setNcg(PetscInt n) { nCG = n; } 
	void setNcgAdj(PetscInt n) { nCGAdj = n; }

	// For backmapping
	bool conv = TRUE;

	std::fstream fpLog;
	PetscViewer viewer;

protected:
    	gmx_int64_t freqUpdate,
    	timeStep,
	Delta, // CG timestep in ps
    	dim,
	cgDim,
	nAtoms,
	nLocalAtoms, // number of atoms per subsystem. TODO: change the identifier to "sub" versus "local"
	nCG,
	nCGx,
	nCGy,
	nCGz,
	nCGAdj,
	neighX, // # of neighbor cells used for doing NNS along the x-direction
	neighY, // # of neighbor cells used for doing NNS along the y-direction
	neighZ, // # of neighbor cells used for doing NNS along the z-direction
	threshold, // minimum number 
	assembleFreq;

	int order, numSS, numSSglob;
	double scale;

    	PetscErrorCode ierr;
	PetscScalar resolution;
    	MPI_Comm comm;

    	PetscInt mpiRank,
	     	 mpiSize;

    	std::vector<PetscScalar> centerOfMass;

    	PetscErrorCode clean();

    	std::string cgMethod;
    	std::vector<Mat> mesoMicroVec,
    			 microMesoVec,
			 kernelVec;

	Mat mesoMicroMap,
            microMesoMap,
            kernel;

	std::map<std::string, ptrMap> fineGrainHash,
				      coarseGrainHash;

	std::map<std::string, ptrMapVelo> coarseGrainVeloHash;

	std::map<std::string, PetscErrorCode (*)(DmsBase&)> initializeHash,
							    updateRefHash; 

	std::fstream& beginLog();

	std::vector< std::list<PetscInt> > fvList;
	std::vector< std::vector<PetscScalar> > grid;
	std::vector< std::vector< std::vector<PetscInt> > > gridID;

private:
    bool debug;
    PetscInt nHistory;
};

#endif

#ifdef __cplusplus
extern "C" {
#endif

void delDmsBase(dmsBasePtr);

// Each public method takes an opaque reference to the object
// that was returned from the above constructor plus the methods parameters.
int dmsCGStep(dmsBasePtr, gmx_int64_t);
bool checkconverge(dmsBasePtr swm);
int constructDmsCoords(dmsBasePtr swm);
int constructDmsVelo(dmsBasePtr swm, const int dmsStep);
dmsBasePtr newDmsBase(const t_state* state, const t_mdatoms* mdatoms,
                      const gmx_mtop_t* top, const t_inputrec* ir, gmx_int64_t dim,
                      gmx_int64_t, int, gmx_int64_t freq, const real dt, const gmx_int64_t,
                      MPI_Comm comm, const int, const float, const int, const PetscInt, const int, char 
		      const*, char* readref, char*, char*);

gmx_bool dmsInitialize(int argc, char* argv[]);

#ifdef __cplusplus
}
#endif

#endif
