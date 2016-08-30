#include "MS_SparseSquare.h"
#include "Utils/Logger.h"

// Construct a 3 * numRows x 3 * numRows sparse square matrix.
MS_SparseSquare::MS_SparseSquare(int numParticles) {
	rows = 3 * numParticles;
	for (int i = 0; i < numParticles; i++) {
		diagonalBlocks.push_back(Matrix3x3::Zero());
	}
}

// Adds the given 3x3 block to this sparse matrix.
void MS_SparseSquare::addBlock(int p1, int p2, Matrix3x3 block) {
	// If it is a diagonal block, add it to whatever is already there.
	if (p1 == p2) {
		diagonalBlocks[p1] += block;
	}
	// Otherwise we know that this block will never be edited again. So
	// add it to the list of non-diagonal blocks.
	else {
		IndexedBlock ib;
		ib.blockRow = p1;
		ib.blockCol = p2;
		ib.block = block;
		offDiagonalBlocks.push_back(ib);
	}
}

// Multiplies this sparse matrix by the given scalar.
void MS_SparseSquare::operator *= (double scalar) {
	for (auto &indexedBlock : offDiagonalBlocks) {
		indexedBlock.block *= scalar;
	}
	for (auto &block : diagonalBlocks) {
		block *= scalar;
	}
}

// Has the effect of adding the mass matrix M defined by the mass vector.
void MS_SparseSquare::addMassAlongDiagonal(dVector &masses) {
	for (int i = 0; i < masses.size(); i++) {
		diagonalBlocks[i] += (Matrix3x3::Identity() * masses[i]);
	}
}

// Solves the system Ax = b, where this sparse matrix is A, and b is the given argument.
// Returns the solution for x.
dVector MS_SparseSquare::solve(dVector b) {
	// First we will have to break the blocks into triplets.
	vector<MTriplet> triplets;

	for (int i = 0; i < rows / 3; i++) {
		int start = i * 3;
		addSparseMatrixDenseBlockToTriplet(triplets, start, start, diagonalBlocks[i]);
	}

	for (auto &indexedBlock : offDiagonalBlocks) {
		int startRow = indexedBlock.blockRow * 3;
		int startCol = indexedBlock.blockCol * 3;
		addSparseMatrixDenseBlockToTriplet(triplets, startRow, startCol, indexedBlock.block);
	}

	// Create a sparse square matrix from our list of triplets
	SparseMatrix A(rows, rows);
	A.setFromTriplets(triplets.begin(), triplets.end());

	// Solve the system with Eigen.
	Eigen::SimplicialLDLT<SparseMatrix, Eigen::Lower> solver;
	solver.compute(A);
	return solver.solve(b);
}