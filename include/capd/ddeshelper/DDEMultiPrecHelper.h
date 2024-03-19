/*
 * DDEMultiPrecHelper.h
 *
 * This file constitutes part of DDEs rigorous integration framework developed
 * in PhD Thesis under supervision of prof. Piotr Zgliczynski:
 *
 * 		"Rigorous Integration of Delay Differential Equations", Jagiellonian University, 2015
 *
 * When using in scientific work please consider citing my articles (preffered),
 * PhD thesis and/or my webpage. For the publications describing the source code
 * please refer to http://scirsc.org/p/papers. Most notable paper up to date is:
 *
 *     Szczelina, R.; Zgliczyński, P.; "Algorithm for rigorous integration of
 *     Delay Differential Equations and the computer-assisted proof of periodic
 *     orbits in the Mackey-Glass equation", http://dx.doi.org/10.1007/s10208-017-9369-5
 *     Foundations of Computational Mathematics (2018), Vol. 18, Iss 6, Pages 1299--1332
 *
 * This work would not be possible without aid and expertise of people involved in
 * CAPD developing library (Computer Assisted Proofs in Dynamics).
 * Please refer to http://capd.ii.uj.edu.pl and consider citing also this library
 * when using those codes in any scientific work.
 *
 * Author: Robert Szczelina, PhD
 * Faculty of Mathematics and Computer Science, Jagiellonian University AND
 * (former) Małopolska Center of Biotechnology, Jagiellonian University
 * email: 	robert.szczelina@uj.edu.pl
 * www: 	scirsc.org
 *
 * This source code is provided under GNU GPL license
 * (v.2 or whatever compatible with CAPD license)
 */

///////////////////////////////////////////////////////////////////////////////////////
// Helper functions to deal with Multiprecision and conversion to/from standard CAPD //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef _DDETOOL_MPIHELPER_H_
#define _DDETOOL_MPIHELPER_H_

#include <capd/vectalg/iobject.hpp>				// helper algorithms for any interval objects (IVectors, intervals, IMatrices, etc...

namespace capd {
namespace ddeshelper {

using namespace capd;
using namespace vectalg;
using namespace intervals;

/** TODO: (NOT URGENT) docs */
interval mpi_to_interval(const MpInterval& a){
	return interval(capd::multiPrec::toDouble(a.leftBound(), MpReal::RoundDown), capd::multiPrec::toDouble(a.rightBound(), MpReal::RoundUp));
};
/** TODO: (NOT URGENT) docs */
MpInterval interval_to_mpi(const interval& a){
	return MpInterval(a.leftBound(), a.rightBound());
};
/** TODO: (NOT URGENT) docs */
template<typename MatrixSpec>
void to_mpi_matrix(const MatrixSpec& inA, MpIMatrix& outA, int precission = 512){
	int old_precision = MpReal::getDefaultPrecision();
	MpReal::setDefaultPrecision(precission);
	int size = inA.numberOfColumns(); // assume square matrix
	outA = MpIMatrix(size, size);
	for (int i = 0; i < size; ++i)
		for (int j = 0; j < size; ++j)
			outA[i][j] = MpInterval(inA[i][j].leftBound(), inA[i][j].rightBound());
	MpReal::setDefaultPrecision(old_precision);
}
/** TODO: (NOT URGENT) docs */
template<typename MatrixSpec>
void from_mpi_matrix(const MpIMatrix& inA, MatrixSpec& outA){
	int size = inA.numberOfColumns(); // assume square matrix
	outA = MatrixSpec(size, size);
	for (int i = 0; i < size; ++i)
		for (int j = 0; j < size; ++j)
			outA[i][j] = mpi_to_interval(inA[i][j]);
}
/** TODO: (NOT URGENT) docs */
template<typename VectorSpec>
void to_mpi_vector(const VectorSpec& inx, MpIVector& outx, int precission = 512){
	int old_precision = MpReal::getDefaultPrecision();
	MpReal::setDefaultPrecision(precission);
	int size = inx.dimension();
	outx = MpIVector(size);
	for (int i = 0; i < size; ++i)
		outx[i] = MpInterval(inx[i].leftBound(), inx[i].rightBound());
	MpReal::setDefaultPrecision(old_precision);
}
/** TODO: (NOT URGENT) docs */
template<typename VectorSpec>
void from_mpi_vector(const MpIVector& inx, VectorSpec& outx){
	int size = inx.dimension();
	outx = VectorSpec(size);
	for (int i = 0; i < size; ++i)
		outx[i] = mpi_to_interval(inx[i]);
}
/** TODO: (NOT URGENT) docs */
template<typename MatrixSpec>
MpIMatrix mpi_inverse(MatrixSpec& A, int precission = 512){
	int old_precision = MpReal::getDefaultPrecision();
	MpReal::setDefaultPrecision(precission);
	MpIMatrix mpA = double_to_mpi_matrix(A, precission);
	mpA.Transpose();
	MpIMatrix mpAinv = capd::matrixAlgorithms::inverseMatrix(mpA);
	mpAinv.Transpose();
	MpReal::setDefaultPrecision(old_precision);
	return mpAinv;
}

// @see: this is to prevent errors when mistakenly mpi_multiply() two MpIMatrices...
void mpi_multiply(MpIMatrix& inA, MpIMatrix& inB, MpIMatrix& outAB, int precission = 512){
	outAB = inA * inB;
}

// TODO: (NOT URGENT) only square matrices...
/** TODO: (NOT URGENT) docs */
template<typename MatrixSpec>
void mpi_multiply(MatrixSpec& inA, MpIMatrix& inB, MatrixSpec& outAB, int precission = 512){
	int old_precision = MpReal::getDefaultPrecision();
	MpReal::setDefaultPrecision(precission);
	int size = inA.numberOfColumns(); // assume square matrix
	MpIMatrix mpA(size, size);
	to_mpi_matrix(inA, mpA, precission);
	MpIMatrix mpAB = mpA * inB;
	outAB = MatrixSpec(size, size);
	from_mpi_matrix(mpAB, outAB);
	MpReal::setDefaultPrecision(old_precision);
}

// TODO: (NOT URGENT) only square matrices...
/** TODO: (NOT URGENT) docs */
template<typename MatrixSpec>
void mpi_multiply(MpIMatrix& inA, MatrixSpec& inB, MatrixSpec& outAB, int precission = 512){
	int old_precision = MpReal::getDefaultPrecision();
	MpReal::setDefaultPrecision(precission);
	int size = inA.numberOfColumns(); // assume square matrix
	MpIMatrix mpB(size, size);
	to_mpi_matrix(inB, mpB, precission);
	MpIMatrix mpAB = inA * mpB;
	outAB = MatrixSpec(size, size);
	from_mpi_matrix(mpAB, outAB);
	MpReal::setDefaultPrecision(old_precision);
}

// TODO: (NOT URGENT) only square matrices...
/** TODO: (NOT URGENT) docs */
template<typename VectorSpec>
void mpi_multiply_Ax(MpIMatrix& inA, VectorSpec& inx, VectorSpec& outAx, int precission = 512){
	int old_precision = MpReal::getDefaultPrecision();
	MpReal::setDefaultPrecision(precission);
	MpIVector mpx(inx.dimension());
	to_mpi_vector(inx, mpx, precission);
	MpIVector mpAx = inA * mpx;
	outAx = VectorSpec(inx.dimension());
	from_mpi_vector(mpAx, outAx);
	MpReal::setDefaultPrecision(old_precision);
}

} // namespace capd
} // namespace ddeshelper

#endif // _DDETOOL_MPIHELPER_H_
