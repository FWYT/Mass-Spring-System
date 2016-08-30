#pragma once

#include "MathLib/P3D.h"
#include "MathLib/V3D.h"
#include "MathLib/Matrix.h"
#include "MS_SparseSquare.h"

class ZeroLengthSpring {
private:
	// Non-stationary particle endpoint of the spring.
	int particle;
	// Spring stiffness.
	double springK;
	// Stationary world position endpoint of the spring.
	P3D restPosition;

	V3D springForce;
	Matrix3x3 dFdxBlock;

public:
	static double globalK;
	static bool useGlobalK;
	ZeroLengthSpring(int p, P3D restPos, double stiffness, const dVector &X);
	int getParticle();
	P3D getRestPosition();
	double getStiffness();

	P3D getCoordinates(int pNum, const dVector& X);

	void computeForces(const dVector& x);
	void addForcesToVector(const dVector& x, dVector& F);

	void compute_dFdx(const dVector& x);
	void addBlocksToMatrix(const dVector& x, MS_SparseSquare &forceJacobian);
};