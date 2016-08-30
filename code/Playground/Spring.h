#pragma once

#include "MathLib/P3D.h"
#include "MathLib/V3D.h"
#include "MathLib/Matrix.h"
#include "MS_SparseSquare.h"

/*
	A class representing a (undamped) spring joining two particles in the system.
	The rest length is considered to be the distance between the two endpoints
	at the time of the spring's creation.
*/

class Spring {
private:
	// Endpoint particles of the spring.
	int particle0, particle1;
	// Rest length of the spring. Computed when the spring is created.
	double restLength;
	// Spring stiffness.
	double springK;

	V3D forces[2];
	Matrix3x3 jacobianBlocks[2][2];

public:
	static double globalK;
	static bool useGlobalK;
	double getStiffness();

	Spring(int p1, int p2, double stiffness, const dVector& X);

	P3D getCoordinates(int pNum, const dVector& X);
	double computeRestLength(const dVector& X);

	void computeForces(const dVector& x);
	void addForcesToVector(const dVector& x, dVector& F);

	void compute_dFdx(const dVector& x);
	void addBlocksToMatrix(const dVector& x, MS_SparseSquare &forceJacobian);

	int p0();
	int p1();
};