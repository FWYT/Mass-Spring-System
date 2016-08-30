#include "Spring.h"
#include "Utils/Logger.h"
#include <math.h>

Spring::Spring(int p0, int p1, double stiffness, const dVector& X) {
	particle0 = p0;
	particle1 = p1;
	restLength = computeRestLength(X);
	springK = stiffness;
}

double Spring::globalK = 1;
bool Spring::useGlobalK = false;

double Spring::getStiffness() {
	if (useGlobalK) {
		return globalK;
	}
	else return springK;
}

// Get the coordinates of a particle.
P3D Spring::getCoordinates(int pNum, const dVector& X) {
	int base = pNum * 3;
	return P3D(X[base], X[base + 1], X[base + 2]);
}

// Compute the rest length of this spring from initial position vector X.
double Spring::computeRestLength(const dVector& X) {
	P3D p1 = getCoordinates(particle0, X);
	P3D p2 = getCoordinates(particle1, X);
	V3D V1(p1, p2);
	//now compute the length of the spring...
	return V1.length();
}

// Compute the forces applied by this spring to its two particles,
// given the current positions of all particles as stored in x.
// forces[0] should store the force applied to particle 0, while
// forces[1] should store the force applies to particle 1.
void Spring::computeForces(const dVector& x) {

	// TODO (Part 2): Compute the forces on both particles.
	forces[0].zero();
	forces[1].zero();

	//particle 0
	int idx0 = p0();
	P3D coor0 = getCoordinates(idx0, x);

	//particle 1 
	int idx1 = p1();
	P3D coor1 = getCoordinates(idx1, x);

	//rest length of spring 
	double restLength = computeRestLength(x);

	//magnitude of vectors |x_i-x_j| 
	//|<u,v>| = sqrt(u^2 + v^2)
	P3D diff = P3D(coor0[0] - coor1[0], coor0[1] - coor1[1], coor0[2] - coor1[2]);
	double mag = sqrt(pow(diff[0], 2) + pow(diff[1], 2) + pow(diff[2], 2));

	//spring stiffness
	double stiffness = getStiffness();
	
	//force for nonzero length spring F = -k((|x_i - x_j|/L) - 1) ((x_i - x_j)/ |x_i - x_j|)
	double scalar = ((-1)*stiffness*((mag / restLength) - 1)) / mag;
	V3D force = V3D(diff[0] * scalar, diff[1] * scalar, diff[2] * scalar);

	//set force vectors
	forces[0] = force;
	forces[1] = V3D(force[0] * -1, force[1] * -1, force[2] * -1);

}

// Add the two forces applied to the endpoints to the force vector.
void Spring::addForcesToVector(const dVector &x, dVector &forceVector) {
	computeForces(x);
	int base0 = 3 * particle0;
	forceVector[base0    ] += forces[0].at(0);
	forceVector[base0 + 1] += forces[0].at(1);
	forceVector[base0 + 2] += forces[0].at(2);

	int base1 = 3 * particle1;
	forceVector[base1    ] += forces[1].at(0);
	forceVector[base1 + 1] += forces[1].at(1);
	forceVector[base1 + 2] += forces[1].at(2);
}

// Compute the four blocks of the Jacobian matrix for this spring's force,
// given the current positions of all particles as stored in x.
// jacobianBlocks[i][j] should be the derivative of force on particle i
// with respect to a change in particle j's position.
void Spring::compute_dFdx(const dVector& x) {

	// TODO (Part 3): Compute the Jacobian matrix blocks for this spring.
	jacobianBlocks[0][0] = Matrix3x3::Zero();
	jacobianBlocks[1][0] = Matrix3x3::Zero();
	jacobianBlocks[0][1] = Matrix3x3::Zero();
	jacobianBlocks[1][1] = Matrix3x3::Zero();
	
	//particle 0
	int idx0 = p0();
	P3D coor0 = getCoordinates(idx0, x);

	//particle 1 
	int idx1 = p1();
	P3D coor1 = getCoordinates(idx1, x);

	//rest length of spring 
	double restLength = computeRestLength(x);

	//magnitude of vectors |x_i-x_j| 
	//|<u,v>| = sqrt(u^2 + v^2)
	V3D diff = coor0 - coor1;
	double mag = sqrt(diff.transpose()*diff);


	//only blocks (1,0) and (0,1) are non zero matrices
	Matrix3x3 uuT = diff.outerProductWith(diff); //(3x1)*(1x3) = (3x3)
	double uTu = diff.transpose().dot(diff); //(1x3)*(3x1) = (1x1) scalar

	double epsilon = (mag / restLength) - 1;

	Matrix3x3 I = Matrix3x3::Zero();
	I(0, 0) = 1;
	I(1, 1) = 1;
	I(2, 2) = 1;

	jacobianBlocks[0][1] = (1 / restLength)*(uuT / uTu) + (I - (epsilon / mag)*(uuT / uTu));
	jacobianBlocks[0][1] *= (-1)*getStiffness();

	//times negative 1 for the other block.
	jacobianBlocks[1][0] = jacobianBlocks[0][1] * (-1);
	
}

// Adds the four 3x3 blocks contributed by this spring to the force Jacobian matrix.
void Spring::addBlocksToMatrix(const dVector &x, MS_SparseSquare &forceJacobian) {
	compute_dFdx(x);
	forceJacobian.addBlock(particle0, particle0, jacobianBlocks[0][0]);
	forceJacobian.addBlock(particle0, particle1, jacobianBlocks[0][1]);
	forceJacobian.addBlock(particle1, particle0, jacobianBlocks[1][0]);
	forceJacobian.addBlock(particle1, particle1, jacobianBlocks[1][1]);
}

int Spring::p0() {
	return particle0;
}

int Spring::p1() {
	return particle1;
}