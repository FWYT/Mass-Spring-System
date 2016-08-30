#pragma once

#include "MathLib/Matrix.h"
#include <vector>

using namespace std;

struct IndexedBlock {
	int blockRow;
	int blockCol;
	Matrix3x3 block;
};

class MS_SparseSquare {
private:
	int rows;
	// The list of blocks off the diagonal stored in this matrix.
	vector<IndexedBlock> offDiagonalBlocks;
	// The blocks on the diagonal. We store these separately because
	// we will need to add blocks together here.
	vector<Matrix3x3> diagonalBlocks;

public:
	MS_SparseSquare(int numParticles);
	void operator *= (double scalar);
	void addMassAlongDiagonal(dVector &masses);
	void addBlock(int p1, int p2, Matrix3x3 block);
	dVector solve(dVector b);
};
