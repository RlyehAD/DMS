#ifndef DMS_HEADERS_H
#define DMS_HEADERS_H

#include "types/enums.h" // must include this first for state.h
#include "types/state.h"
#include "types/atoms.h"
#include "types/topology.h"
#include "types/simple.h"
#include "gmx_fatal.h"
#include "mdatoms.h"

#ifdef __cplusplus

// PETSc
#include <petsc.h>
#include <petscksp.h>
#include <petscmat.h>
#include <petscsys.h>
#include <petsclog.h>

// C++ stl
#include <iostream>
#include <vector>
#include <memory>
#include <list>
#include <functional>
#include <algorithm>
#include <map>
#include <stdexcept>
#include <fstream>
#include <string>
#include <sstream>
#include <iterator>

/* the idea behind this inline function is to check on errors in functions that
 * cannot return anything such as constructors and destructors.
 */
inline PetscErrorCode DMS_CHKERRQ(PetscErrorCode ierr) { CHKERRQ(ierr); }
PetscScalar CalcMHWScore(std::vector<PetscScalar>);
typedef std::vector<Vec> CVec;

//gmx_fatal(FARGS, "Error from DMS %d", 1);

#endif
#endif
