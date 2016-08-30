#pragma once

#include "Spring.h"
#include "ZeroLengthSpring.h"
#include "MS_SparseSquare.h"
#include <vector>
#include <list>
#include "GUILib/GLMesh.h"

using namespace std;

struct Particle {
	P3D position;
	P3D velocity;
	double mass;
};

using namespace Eigen;

class ParticleSystem {
private:
	// Vectors containing particle states for the system.
	dVector positions;
	dVector velocities;
	dVector masses;
	int numParticles;

	// A list of all springs in the system.
	vector<Spring> springs;

	// A list of all zero length springs in the system.
	list<ZeroLengthSpring> zeroLengthSprings;

	// A list of flags for each particle specifying if it's currently pinned.
	vector<bool> particlePinned;

	// Vectors to pass to OpenGL for drawing.
	// Each time step, the relevant data are copied into these lists.
	vector<double> positionArray;
	vector<unsigned int> pointsIndexArray;
	vector<unsigned int> edgesIndexArray;
	vector<double> zlSpringPositionArray;
	int count;

	GLMesh* mesh;
	
public:
	ParticleSystem(vector<Particle>& particles);
	void addSpring(int p1, int p2, double stiffness);
	void addZeroLengthSpring(int p, P3D x0, double stiffness);
	void deleteZeroLengthSpring(int p);
	void togglePin(int p);
	P3D getPositionOf(int i);
	int particleCount();
	dVector computeForceVector(const dVector &x, const dVector &v);
	MS_SparseSquare computeForceJacobian(const dVector &x, const dVector &v);
	dVector computeAccelerations(const dVector &x, const dVector &v, const dVector &m);
	void integrate_FE(double delta);
	void integrate_SE(double delta);
	dVector newtonStep(double delta, dVector &positions, dVector &velocities, dVector &masses, dVector &v_guess);
	void integrate_BE(double delta);

	// Functions for display and interactivity
	void drawParticleSystem();
	void setPosition(int i, P3D x);
	void setVelocity(int i, V3D v);
	void setMesh(GLMesh* m);

	// Whether or not we should draw springs and particles as lines and dots respectively.
	static bool drawSprings;
	static bool drawParticles;
	static bool drawZeroLengthSprings;
	static bool drawFaces;
	static bool enableGravity;
};
