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

#ifndef _CAPD_DDES_DDESOLUTIONCURVE_DOTIMPL_HPP_
#define _CAPD_DDES_DDESOLUTIONCURVE_DOTIMPL_HPP_

#include <capd/ddes/DDESolutionCurve.h>
#include <capd/ddes/DDECommon.hpp>
#include <capd/ddes/DDEForwardTaylorCurvePiece.hpp>

////////////////////////////////////////////////////////////////////////
// Below are implementations of more involved functions               //
// (to be separated as usual from the declarations in the .h file     //
// For more complete documentation/explanation please consult .h file //
////////////////////////////////////////////////////////////////////////

namespace capd{
namespace ddes{

template<typename SetSpec>
typename DDESolutionCurve<SetSpec>::ScalarType
DDESolutionCurve<SetSpec>::dot(DDESolutionCurve<SetSpec>::VectorType const& v) const{
	// TODO: (NOT URGENT) dim check?
	size_type vd = v.dimension();
	size_type N0 = storageN0();
	size_type d = dimension();
	MatrixType CT(N0, vd);
	MatrixType BT(vd, vd);
	VectorType r(vd);
	VectorType x(vd);
	MatrixType cC = m_valueAtCurrent.get_C();
	MatrixType cB = m_valueAtCurrent.get_B();
	VectorType cr = m_valueAtCurrent.get_r();
	VectorType cx = m_valueAtCurrent.get_x();
	size_type I = 0;
	size_type J = 0;
	for (size_type i = 0; i < d && I < vd; ++i, ++I){
		x[J + i] = cx[i];
		r[J + i] = cr[i];
		for (size_type j = 0; j < N0; ++j)
			CT[j][J + i] = cC[i][j];
		for (size_type j = 0; j < d; ++j)
			BT[J + j][J + i] = cB[i][j];
	}
	J += d;
	for (auto ijet = rbegin(); ijet != rend() && I < vd; ++ijet){
		for (auto icoeff = (*ijet)->beginJet(); icoeff != (*ijet)->endJet() && I < vd; ++icoeff){
			cC = icoeff->get_C();
			cB = icoeff->get_B();
			cr = icoeff->get_r();
			cx = icoeff->get_x();
			for (size_type i = 0; i < d && I < vd; ++i, ++I){
				x[J + i] = cx[i];
				r[J + i] = cr[i];
				for (size_type j = 0; j < N0; ++j)
					CT[j][J + i] = cC[i][j];
				for (size_type j = 0; j < d; ++j)
					BT[J + j][J + i] = cB[i][j];
			}
			J += d;
		}
	}
	return v * x + (CT * v) * (*m_r0) + (BT * v) * r;
}

template<typename SetSpec>
template<typename JetSection>
typename DDESolutionCurve<SetSpec>::ScalarType
DDESolutionCurve<SetSpec>::dot(JetSection const& v) const{
	// TODO: (NOT URGENT) dim check?
	size_type vd = v.storageDimension();
	size_type N0 = storageN0();
	size_type d = dimension();
	MatrixType CT(N0, vd);
	MatrixType BT(vd, vd);
	VectorType r(vd);
	VectorType x(vd);
	MatrixType cC = m_valueAtCurrent.get_C();
	MatrixType cB = m_valueAtCurrent.get_B();
	VectorType cr = m_valueAtCurrent.get_r();
	VectorType cx = m_valueAtCurrent.get_x();
	size_type I = 0;
	size_type J = 0;
	for (size_type i = 0; i < d && I < vd; ++i, ++I){
		x[J + i] = cx[i];
		r[J + i] = cr[i];
		for (size_type j = 0; j < N0; ++j)
			CT[j][J + i] = cC[i][j];
		for (size_type j = 0; j < d; ++j)
			BT[J + j][J + i] = cB[i][j];
	}
	J += d;
	auto vjet = v.rbegin();
	auto cjet = rbegin();
	for (; cjet != rend() && vjet != v.rend(); ++cjet, ++vjet){
		size_type order = vjet->order();
		auto icoeff = (*cjet)->beginJet();
		for (size_type k = 0; k <= order && I < vd; ++icoeff, ++k){
			if (icoeff == (*cjet)->endJet())
				throw std::logic_error("DDESolutionCurve::dot(JetSection): curve has not enough order at some jet to match the section.");
			cC = icoeff->get_C();
			cB = icoeff->get_B();
			cr = icoeff->get_r();
			cx = icoeff->get_x();
			for (size_type i = 0; i < d && I < vd; ++i, ++I){
				x[J + i] = cx[i];
				r[J + i] = cr[i];
				for (size_type j = 0; j < N0; ++j)
					CT[j][J + i] = cC[i][j];
				for (size_type j = 0; j < d; ++j)
					BT[J + j][J + i] = cB[i][j];
			}
			J += d;
		}
	}
	VectorType vv(v);
	return vv * x + (CT * vv) * (*m_r0) + (BT * vv) * r;
}


} // namespace ddes
} // namespace capd

#endif /* _CAPD_DDES_DDESOLUTIONCURVE_DOTIMPL_HPP_ */
