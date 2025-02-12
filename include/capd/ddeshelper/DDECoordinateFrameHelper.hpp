/*
 * >>PUT YOUR FILENAME HERE<<.hpp
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

#ifndef _CAPD_DDECOORDINATEFRAMEHELPER_HPP_
#define _CAPD_DDECOORDINATEFRAMEHELPER_HPP_

#include <capd/ddeshelper/DDECoordinateFrameHelper.h>


namespace capd {
namespace ddeshelper {

/**
 * Implements finding a coordinate frame for a dynamical system/stochastic process
 * with the Karhunen–Loève Transform / PCA (see Kosambi-Karhunen–Loève theorem).
 *
 * The standard $\R^M$ scalar product is used.
 *
 * @param M is the dimension of the vectors in U (TODO: (NOT URGENT< RETHINK) - remove it as unnecessary?)
 * @param U is the set of vectors from the process (time snapshots).
 * @param avg (out) here will be the avg pasition of U
 * @param coords (out) here will be the coordinates around avg for the process (the eigenfunction of the correlation matrix)
 * @param eigen (out) here will be the eigenvalues associated to the respective coordinates in the matrix (dominant first)
 * @param eigen normalize (optional, default: true), if set to true, then the columns of coords will be normalized
 */
template<typename MatrixSpec, typename VectorSpec>
void correlationCoords(
		int M,
		std::vector<VectorSpec> U,
		VectorSpec& avg,
		MatrixSpec& coords,
		VectorSpec& eigen,
		bool normalize=true){
	typedef MatrixSpec MatrixType;
	typedef VectorSpec VectorType;

	avg = VectorType(M);
	for (auto u = U.begin(); u != U.end(); ++u)
		avg += (*u);
	avg /= (double)U.size();

	for (auto u = U.begin(); u != U.end(); ++u)
		(*u) -= avg;

	MatrixType A(M, M);
	for (auto u = U.begin(); u != U.end(); ++u)
		for (int i = 0; i < M; ++i)
			for (int j = 0; j < M; ++j)
				A[i][j] += (*u)[i] * (*u)[j];
	A /= U.size();

	VectorType valIm(M);
	MatrixType vecIm(M, M);
	eigen = VectorType(M);
	coords = MatrixType(M, M);
	capd::alglib::computeEigenvaluesAndEigenvectors(A, eigen, valIm, coords, vecIm);

	bool good = true;
	std::ostringstream info;
	info << "correlationCoords(): eigenvalue with nonzero imagine part\n";
	for (int i = 0; i < M; ++i){
		if (std::abs(valIm[i]) > 1e-15){
			info << "\t for i = " << i << " " << eigen[i] << " + " << valIm[i] << "i\n";
			good = false;
		}
	}
	if (!good) throw std::logic_error(info.str());

	if (normalize)
		for (int i = 0; i < M; ++i)
			coords.column(i).normalize();
}

} // namespace ddeshelper
} // namespace capd

#endif /* _CAPD_DDECOORDINATEFRAMEHELPER_HPP_ */
