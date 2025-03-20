/*
 * ddeshelperlib.h
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

#ifndef _CAPD_DDES_DDESHELPERSLIB_H_
#define _CAPD_DDES_DDESHELPERSLIB_H_

//#ifdef _CAPD_MPCAPD_H_
//#include <capd/ddeshelper/DDEMultiPrecHelper.h>
//#endif

#include <capd/ddeshelper/DDEHelperCommon.h>
#include <capd/ddeshelper/DDESystemOSHelper.h>
#include <capd/ddeshelper/DDEHelperDrawing.hpp>
#include <capd/ddeshelper/DDECompareHelper.h>
#include <capd/ddeshelper/DDECoordinateFrameHelper.hpp>
// user should decide by himself
//#include <capd/ddeshelper/DDEHelperNonrigorous.hpp>
//#include <capd/ddeshelper/DDEHelperRigorous.hpp>


// TODO: (NOT IMPORTANT): Move to some other file.
namespace capd{
namespace ddeshelper{

/**
 * Works on two solutions (collections of jets + value at t0) of the same kind.
 * Makes one solution into another by dumping higher order terms in jets
 *
 * ***BEWARE***: (1) it copies the value at t0, (2) the output jet controls
 * how many jets are copied. If not enough jets in v then ERROR
 * (seg. fault most probably)!
 */
template<typename SolType>
void copyReduce(const SolType &v, SolType &u) {
	u.setValueAtCurrent(v.getValueAtCurrent());
	auto jetu = u.begin();
	auto jetv = v.begin();
	for (; jetu != u.end(); ++jetu, ++jetv){
		auto coeffu = (*jetu)->begin();
		auto coeffv = (*jetv)->begin();
		for (; coeffu != (*jetu)->end(); ++coeffu, ++coeffv){
			(*coeffu) = (*coeffv);
		}
	}
}

} // namespace ddeshelper
} // namespace capd

#endif /* _CAPD_DDES_DDESHELPERSLIB_H_ */
