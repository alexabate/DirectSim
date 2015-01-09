/**
 * @file  matrix.h
 * @brief Matrix routines, SVD
 *
 * Could add more information here I think
 *
 * @author Alex Abate
 * Contact: abate@email.arizona.edu
 *
 * Created on: 2012
 * @date 2012
 *
 */
#ifndef  IntfLapack_H_SEEN
#ifndef MATRIX_SEEN
#define MATRIX_SEEN

#include <math.h>
#include <vector>

// sophya
#include "sopnamsp.h"
#include "machdefs.h"
#include "tvector.h"
#include "cspline.h"

namespace SOPHYA {


struct SVD 
	{
	/** Object for singular value decomposition of a matrix A, and related functions */
	int m,n;
	TMatrix<double> u,v;
	TVector<double> w;
	double eps,tsh;

	/** Constructor */
	SVD(TMatrix<double> &a); 

	/** Calculate the inverse */
	TMatrix<double> matrix_inv_svd();

	/** Solve with (apply the pseudoinverse to) one or more right-hand sides */
	void solve(TVector<double> &b, TVector<double> &x, double thresh);
	
	/** Solve with (apply the pseudoinverse to) one or more right-hand sides */
	void solve(TMatrix<double> &b, TMatrix<double> &x, double thresh);

	// Quantities associated with the range and nullspace of A
//	int rank(double thresh);
//	int nullity(double thresh);
//	TMatrix<double> range(double thresh);
//	TMatrix<double> nullspace(double thresh);

	// Return reciprocal of the condition number of A
//	double inv_condition()
//		{ return (w(0) <= 0. || w(n-1) <= 0.) ? 0. : w(n-1)/w(0); }

	void decompose();
	void reorder();
	double pythag(const double a, const double b);
};

TArray<double> Transpose(TArray<double>);
TMatrix<double> Mult(TMatrix<double>,TMatrix<double>);
TArray<double> Mult(TArray<double>,TArray<double>);

//TArray<double> columnRowMult(vector<double> columnVector,vector<double> rowVector);
//double vectorMagnitude(vector<double> vect);

}  // End namespace SOPHYA

#endif
#endif
