/* aggs.h
 * 
 * PKDGRAV Source Code
 * 
 * Author: Kenneth W. Flynn
 *         flynnk@astro.umd.edu
 * Mods:   Derek C. Richardson
 *         dcr@astro.umd.edu
 *
 * Modified: 01/28/01, DCR: 07/10/02
 * 
 * Aggregate expansion to PKDGRAV.
 *
 * This is the main header file for the aggregate expansion to pkdgrav.  This
 *  file and other associated changes extend the pkdgrav code to handle one or
 *  more aggregate structures.  Each aggregate is a collection of particles
 *  that do not mutually interact but that can interact with outside particles
 *  via collisions and/or gravity torques.  Rotation of the aggregate is
 *  computed using Euler's equations.
 *
 * To enable this functionality, define the AGGS flag.
 *
 * Aggregates make few changes to the original PKD code.  During the initial
 *  setup, the aggregrate structures defined below need to be filled out.
 *  Calculation of gravity is unaffected, we merely make aggregrate particles
 *  "inactive" during a kick, making use of pre-exisiting functionality.
 *  During drift, we update particle positions in the aggregate "by hand"
 *  under the influence of their rotation, computed via the Euler equations.
 *  Collisions are computed by .... *DEBUG* update!
 */
 
#ifndef __AGGS_H
#define __AGGS_H

#ifdef AGGS

#include "pkd.h"
#include "linalg.h"
#include "collision.h"

#define AGGS_INIT_BUFF_SIZE 16 /* allocate this many aggs minimum */

/* Aggregate macro prototypes */
int IS_AGG(PARTICLE *p);
int AGG_ID(PARTICLE *p);
int AGG_SET_ID(PARTICLE *p,int ID);
int COLLIDER_IS_AGG(COLLIDER *c);
int COLLIDER_AGG_ID(COLLIDER *c);
int COLLIDER_AGG_SET_ID(COLLIDER *c,int ID);

/* Particle is an aggregate if its original index is negative */
#define IS_AGG(p) ((p)->iOrgIdx < 0)

/* Aggregate index = -1 - (particle original index) */
#define AGG_ID(p) (-1 - (p)->iOrgIdx)

#define AGG_SET_ID(p,ID) {(p)->iOrgIdx = -1 - (ID);}

/* Same thing for colliders */
#define COLLIDER_IS_AGG(c) ((c)->id.iOrgIdx < 0)
#define COLLIDER_AGG_ID(c) (-1 - (c)->id.iOrgIdx)
#define COLLIDER_AGG_SET_ID(c,ID) {(c)->id.iOrgIdx = -1 - (ID);}

/* Constant in I = 2/5 mr^2 (for spheres).  Change for other bodies. */
#define GAMMA 0.4

/***************************** Structures *****************************/

/* Main structure that defines the aggregate. */
typedef struct 
{
	/* Flag showing whether any particles are assigned to this agg. */
	int bAssigned;

	/* COM position of the aggregate. */
	Vector r_com;

	/* COM velocity of the aggregate. */
	Vector v_com;

	/* COM acceleration of the aggregate. */
	Vector a_com;

	/* Torques acting on the aggregate. */
	Vector torque;
 
	/* Angular velocity of the aggregate. */
	Vector omega;
 
	/*
	 ** Moments of inertia for the aggregate.  Note this is not a spatial
	 ** vector, but a collection of three values.
	 */
	Vector moments;
 
	/*
	 ** Transformation matrix from body to space frame, called lambda. 
	 ** This is simply a matrix whose columns are the principal axes of
	 ** inertia. A call to matrixTranspose() will yield the inverse
	 ** conversion.
	 */
	Matrix lambda;

	/* Time of last update (within step interval) for "drifting" */

	double dLastUpdate;

	} Aggregate;

/***************************** Functions *****************************/

void pkdAggsFind(PKD pkd,int *iMaxIdx);

void pkdAggsConfirm(PKD pkd,int iAggIdx,int *bAssigned);

void pkdAggsMerge(PKD pkd,int iOldID,int iNewID);

void pkdAggsGetCOM(PKD pkd,int iAggIdx,Scalar *m,Vector mr,Vector mv);

void pkdAggsGetAxes(PKD pkd,int iAggIdx,Vector r_com,Vector v_com,
					Matrix I,Vector L);

void pkdAggsToBodyAxes(PKD pkd,int iAggIdx,Matrix spaceToBody);

void pkdAggsSetSpacePos(PKD pkd,int iAggIdx,Vector r_com,Matrix lambda);

void pkdAggsSetSpaceVel(PKD pkd,int iAggIdx,Vector v_com,Vector omega,
						Matrix lambda);

void pkdAggsAccel(PKD pkd,int iAggIdx,Scalar *m,Vector ma);

void pkdAggsTorque(PKD pkd,int iAggIdx,Vector r_com,Vector a_com,
				   Vector torque);

void pkdAggsActivate(PKD pkd);

void pkdAggsDeactivate(PKD pkd);

/** Function that returns the values for the derivatives of the Euler 
 *   equations of motion.  For more information, refer to the accompanying
 *   document.
 *
 *  Parameters (in):
 *   t - Time to evaluate at (ignored)
 *   vars - Current value of variables being integrated.  In detail:
 *            0 = omega[0]
 *            1 = omega[1]
 *            2 = omega[2]
 *            3 = q1[0]
 *            4 = q1[1]
 *            5 = q1[2]
 *            6 = q2[0]
 *            7 = q2[1]
 *            8 = q2[2]
 *            9 = q3[0]
 *           10 = q3[1]
 *           11 = q3[2]
 *   agg_as_void - The aggregate structure, used to obtain the principal 
 *                  moments of inertia
 *
 *  Parameters (out):
 *   derivs - Values for the derivatives of the aforementioned variables
 */
void aggsEulerDerivs(FLOAT t,FLOAT vars[],void* agg_as_void,FLOAT derivs[]);

/** Functions passed to aggsRungeStep() must fit the criterion laid out here.
 *   We provide a description of the parameters here for reference.  Functions
 *   that implement this are essentially differential equations.  The goal
 *   of this function is to identify the value of the differential equations
 *   at the specified point in parameters space.
 *
 *  Parameters (in):
 *   x - Current value of the independent variable
 *   y_vals - Array of values of the dependent variables
 *   user_data - A pointer to the information passed into rk4step, which can
 *                be used to access arbitrary data
 *
 *  Parameters (out):
 *   deriv_vals - The evaluated values for the differential equations
 */
typedef void (*aggsRungeDerivs)(FLOAT x,FLOAT y_vals[],void* user_data,
								FLOAT deriv_vals[]);

/** Take an rk4 step.
 *
 *  Parameters (in):
 *   step_size - How big a step to take
 *   x - Current value of the independent variable
 *   y_vals - Current values for the dependent variables
 *   n - Number of dependent variables
 *   user_data - Arbitrary data submitted by the user; aggsRungeStep() does
 *                not use this data directly; instead it passes it to the
 *                aggsRungeDerivs function
 *   func - Function to evaluate the values of the derivatives at a given
 *           set of values for the independent and dependent variables
 *
 *  Parameters (out):
 *   new_x - New value of the independent variable
 *   new_y_vals - New value for the dependent variables
 */
void aggsRungeStep(FLOAT step_size,FLOAT x,FLOAT* y_vals,int n,void* user_data,
				   aggsRungeDerivs func,FLOAT* new_x,FLOAT* new_y_vals);

#endif /* AGGS */

#endif
