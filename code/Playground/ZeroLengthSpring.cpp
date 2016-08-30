#include "ZeroLengthSpring.h"
#include "Utils/Logger.h"
#include <iostream>

ZeroLengthSpring::ZeroLengthSpring(int p, P3D restPos, double stiffness, const dVector &X) {
	restPosition = restPos;
	springK = stiffness;
	particle = p;
}

double ZeroLengthSpring::globalK = 100;
bool ZeroLengthSpring::useGlobalK = false;

// Get the coordinates of a particle.
P3D ZeroLengthSpring::getCoordinates(int pNum, const dVector& X) {
	int base = pNum * 3;
	return P3D(X[base], X[base + 1], X[base + 2]);
}

double ZeroLengthSpring::getStiffness() {
	if (useGlobalK) {
		return globalK;
	}
	else return springK;
}

int ZeroLengthSpring::getParticle() {
	return particle;
}

P3D ZeroLengthSpring::getRestPosition() {
	return restPosition;
}

// Compute the force applied by this spring to its particle,
// given the current positions of all particles as stored in x.
// The result should be stored in the field springForce.
void ZeroLengthSpring::computeForces(const dVector& x) {

	// TODO (Part 2): Compute the force on the particle and store it in springForce.
	springForce.zero();

	//for zero length spring, F = -k(x_i - x_j), k = stiffness
	double stiffness = getStiffness();
	//stationary endpoint position
	P3D stationary = getRestPosition();
	//non stationary endpoint position
	int idx = getParticle();
	P3D nonStationary = getCoordinates(idx, x);

	//Logger::consolePrint("Particle number: %d", idx);
	//Logger::consolePrint("diff: %d %d %d", (nonStationary[0] - stationary[0]), (nonStationary[1] - stationary[1]), (nonStationary[2] - stationary[2]));
	//Logger::consolePrint("stiffness: %d", stiffness);

	for (int i = 0; i < 3; i++)
	{
		//F = -k(x_i - x_j)
		springForce[i] = (-1)*stiffness*(nonStationary[i] - stationary[i]);	
		if (i == 1)
		{
			//gravity?
			springForce[i] += 9.8;
		}
		//Logger::consolePrint("force[%d] = %d", i, springForce[i]);
	}

	
	

}

// Add the force applied to the endpoint to the force vector.
void ZeroLengthSpring::addForcesToVector(const dVector& x, dVector& F) {
	computeForces(x);
	int base = 3 * particle;
	F[base    ] += springForce.at(0);
	F[base + 1] += springForce.at(1);
	F[base + 2] += springForce.at(2);
}

int countt = 0;

// Compute the single block of the Jacobian matrix for this spring's force,
// given the current positions of all particles as stored in x.
// The result should be stored in the field dFdxBlock.
void ZeroLengthSpring::compute_dFdx(const dVector& x) {

	// TODO (Part 3): Compute the Jacobian matrix for zero-length springs and store it in dFdxblock.
	dFdxBlock = Matrix3x3::Zero();
	

	//df/dx = -kI, I=identity matrix
	dFdxBlock(0, 0) = 1.;
	dFdxBlock(1, 1) = 1;
	dFdxBlock(2, 2) = 1;

	dFdxBlock *= (-1)*getStiffness();
		
}

// Adds the 3x3 block contributed by this spring to the force Jacobian matrix.
void ZeroLengthSpring::addBlocksToMatrix(const dVector& x, MS_SparseSquare &forceJacobian) {
	compute_dFdx(x);
	
	forceJacobian.addBlock(particle, particle, dFdxBlock);
}