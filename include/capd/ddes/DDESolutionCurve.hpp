/*
 * DDESolutionCurve.hpp
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

#ifndef _CAPD_DDES_DDESOLUTIONCURVE_HPP_
#define _CAPD_DDES_DDESOLUTIONCURVE_HPP_

#include <capd/ddes/DDESolutionCurve.h>
#include <capd/ddes/DDECommon.hpp>

////////////////////////////////////////////////////////////////////////
// Below are implementations of more involved functions               //
// (to be separated as usual from the declarations in the .h file     //
// For more complete documentation/explanation please consult .h file //
////////////////////////////////////////////////////////////////////////

#include <capd/ddes/DDESolutionCurve_IO.hpp>
#include <capd/ddes/DDESolutionCurve_PiecesManip.hpp>
#include <capd/ddes/DDESolutionCurve_DoubletonImpl.hpp>
#include <capd/ddes/DDESolutionCurve_Move.hpp>
#include <capd/ddes/DDESolutionCurve_dotImpl.hpp>

namespace capd{
namespace ddes{

template<typename SetSpec>
void hull(DDESolutionCurve<SetSpec> const& X, DDESolutionCurve<SetSpec> const& Y, DDESolutionCurve<SetSpec>& Z){
	if (X.length() != Y.length())
		throw std::logic_error("capd::ddes::hull(): X and Y curves have different lengths.");
	if (X.length() != Z.length())
		throw std::logic_error("capd::ddes::hull(): output curve has different length.");
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

	// we do not care about different time points at Jets.
	auto xit = X.begin();
	auto yit = Y.begin();
	auto zit = Z.begin();
	for ( ; xit != X.end(); ++xit, ++yit, ++zit){
		if ((*xit)->order() != (*yit)->order())
			throw std::logic_error("capd::ddes::hull(): X and Y have different order at some point");
		if ((*xit)->order() != (*zit)->order())
			throw std::logic_error("capd::ddes::hull(): output storage has different order at some point");
		auto xcend = (*xit)->endJet();
		auto xcoeff = (*xit)->beginJet();
		auto ycoeff = (*yit)->beginJet();
		auto zcoeff = (*zit)->beginJet();
		for ( ; xcoeff != xcend; ++xcoeff, ++ycoeff, ++zcoeff)
			capd::ddes::hull(*xcoeff, *ycoeff, *zcoeff);
	}
	capd::ddes::hull(X.getValueAtCurrent(), Y.getValueAtCurrent(), Z.getValueAtCurrent());
}

} // namespace ddes
} // namespace capd

#endif /* _CAPD_DDES_DDESOLUTIONCURVE_HPP_ */
