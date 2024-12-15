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

#ifndef _CAPD_DDES_DDESOLUTIONCURVE_PIECESMANIP_HPP_
#define _CAPD_DDES_DDESOLUTIONCURVE_PIECESMANIP_HPP_

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
DDESolutionCurve<SetSpec>&
DDESolutionCurve<SetSpec>::addPiece(CurvePieceType* newPiece, bool passOwnership){
	if (newPiece->storageN0() != storageN0()){
		std::ostringstream info;
		info << "DDESolutionCurve::addPiece(): the new piece has different N_0.";
		info << "Should be: " << storageN0() << ", is " << newPiece->storageN0() << ".";
		throw std::logic_error(info.str());
	}
	if (!passOwnership)
		newPiece = new CurvePieceType(*newPiece);
	m_pieces.push_back(newPiece);
	m_pieceOwner.push_back(passOwnership);
	newPiece->setT0(currentTime());
	newPiece->set_r0(m_r0);
	++m_t_current;
	return *this;
}
template<typename SetSpec>
DDESolutionCurve<SetSpec>&
DDESolutionCurve<SetSpec>::addPiece(CurvePieceType const & newPiece){
	return this->addPiece(new CurvePieceType(newPiece), true);
}

template<typename SetSpec>
DDESolutionCurve<SetSpec>&
DDESolutionCurve<SetSpec>::addPastPiece(CurvePieceType* newPiece, bool passOwnership){
	if (newPiece->storageN0() != storageN0())
		throw std::logic_error("DDESolutionCurve::addPastPiece(): the new piece has different N_0.");
	if (!passOwnership)
		newPiece = new CurvePieceType(*newPiece);
	m_pieces.push_front(newPiece);
	m_pieceOwner.push_front(passOwnership);
	newPiece->setT0(pastTime() - 1);
	newPiece->set_r0(m_r0);
	return *this;
}

template<typename SetSpec>
DDESolutionCurve<SetSpec>&
DDESolutionCurve<SetSpec>::addPastPiece(const CurvePieceType& newPiece){
	return this->addPastPiece(new CurvePieceType(newPiece), true);
}

template<typename SetSpec>
DDESolutionCurve<SetSpec>&
DDESolutionCurve<SetSpec>::setValueAtCurrent(const SetType& value){
	if (value.storageN0() != storageN0())
		throw std::logic_error("DDESolutionCurve::setValueAtCurrent(set): provided value has incompatible N0 dimension.");
	if (value.dimension() != dimension())
		throw std::logic_error("DDESolutionCurve::setValueAtCurrent(set): provided value has incompatible dimension.");
	this->m_valueAtCurrent = value;
	if (!value.common_r0(m_r0))
		*m_r0 = intervalHull(*m_r0, value.get_r0());
	return *this;
}

template<typename SetSpec>
DDESolutionCurve<SetSpec>&
DDESolutionCurve<SetSpec>::setValueAtCurrent(const VectorType& value){
	if (value.dimension() != dimension())
		throw std::logic_error("DDESolutionCurve::setValueAtCurrent(vector): provided value has incompatible dimension.");
	this->m_valueAtCurrent = SetSpec(value, m_r0); return *this;
}

} // namespace ddes
} // namespace capd

#endif /* _CAPD_DDES_DDESOLUTIONCURVE_PIECESMANIP_HPP_ */
