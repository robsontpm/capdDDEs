/*
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

#ifndef _CAPD_DDES_DOUBLETONINTERFACE_HPP_
#define _CAPD_DDES_DOUBLETONINTERFACE_HPP_

#include <capd/ddes/storage/DoubletonInterface.h>
#include <capd/ddes/DDECommon.hpp>

////////////////////////////////////////////////////////////////////////
// Below are implementations of more involved functions               //
// (to be separated as usual from the declarations in the .h file     //
// For more complete documentation/explanation please consult .h file //
////////////////////////////////////////////////////////////////////////

namespace capd{
namespace ddes{

template<typename SetSpec>
void hull(SetSpec const& X, SetSpec const& Y, SetSpec& Z){
	typedef SetSpec SetType;
	typedef typename SetType::VectorType VectorType;
	typedef typename SetType::MatrixType MatrixType;
	typedef typename SetType::size_type size_type;

	if (X.dimension() != Y.dimension())
		throw std::logic_error("capd::ddes::hull(): X and Y have different dimensions!");
	if (X.dimension() != Z.dimension())
		throw std::logic_error("capd::ddes::hull(): output storage has different dimensions!");
	if (X.storageN0() != Y.storageN0())
		throw std::logic_error("capd::ddes::hull(): X and Y have different N0 dimensions!");
	if (X.storageN0() != Z.storageN0())
		throw std::logic_error("capd::ddes::hull(): output storage has different N0 dimensions!");
	if (X.storageDimension() != Y.storageDimension())
		throw std::logic_error("capd::ddes::hull(): X and Y have different storage dimensions!");
	if (X.storageDimension() != Z.storageDimension())
		throw std::logic_error("capd::ddes::hull(): output storage has different storage dimensions!");

	size_type sd = X.storageDimension();
	size_type N0 = X.storageN0();
	VectorType r0(N0); capd::vectalg::intervalHull(X.get_r0(), Y.get_r0(), r0);
	VectorType x(sd); capd::vectalg::intervalHull(X.get_x() + X.get_B() * X.get_r(), Y.get_x() + Y.get_B() * Y.get_r(), x);
	MatrixType Id(sd, sd); Id.setToIdentity();
	MatrixType C(sd, N0); capd::vectalg::intervalHull(X.get_C(), Y.get_C(), C);
	MatrixType S(sd, N0); capd::vectalg::split(C, S);
	x += S * r0;
	VectorType r(sd); capd::vectalg::split(x, r);
	Z.set_x(x);
	Z.set_Cr0(C, r0);
	Z.set_B(Id);
	Z.set_Binv(Id);
	Z.set_r(r);
}

} // namespace ddes
} // namespace capd

#endif /* _CAPD_DDES_DOUBLETONINTERFACE_HPP_ */
