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

#ifndef _CAPD_DDES_DDECURVEPIECE_HPP_
#define _CAPD_DDES_DDECURVEPIECE_HPP_

#include <capd/ddes/DDECommon.hpp>
#include <capd/ddes/DDEForwardTaylorCurvePiece.h>
#include <capd/ddes/storage/GenericJet.hpp>

////////////////////////////////////////////////////////////////////////
// Below are implementations of more involved functions               //
// (to be separated as usual from the declarations in the .h file     //
// For more complete documentation/explanation please consult .h file //
////////////////////////////////////////////////////////////////////////

namespace capd{
namespace ddes{

template<typename TimePointSpec, typename SetSpec, bool isInterval>
DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>&
DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::operator=(
		DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval> const & orig
)
{
	// we need to handle memory allocations carefuly
	deallocate();
	m_dimension = orig.m_dimension;
	m_order = orig.m_order;
	m_r0_owner = orig.m_r0_owner;
	m_Xi_owner = orig.m_Xi_owner;
	if (orig.m_r0_owner){ allocateR0(orig.storageN0()); *m_r0 = *orig.m_r0; } else { m_r0 = orig.m_r0; }
	if (orig.m_Xi_owner){ allocateXi(); *m_Xi = *orig.m_Xi; } else { m_Xi = orig.m_Xi; }
	allocateJet();
	auto thisJetIt = beginJet();
	auto origJetIt = orig.beginJet();
	for ( ; origJetIt != orig.endJet(); ++origJetIt, ++thisJetIt) *thisJetIt = *origJetIt;
	updateCommonR0();
	return *this;
}

template<typename TimePointSpec, typename SetSpec, bool isInterval>
DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::DDEForwardTaylorCurvePiece(
		TimePointType t0):
	m_t0(t0), m_dimension(0), m_order(0),
	m_jet_at_t0(0),
	m_r0(0), m_r0_owner(true),
	m_Xi(0), m_Xi_owner(true)
{
	allocate();
	setAsConstant(VectorType());
	updateCommonR0();
}

template<typename TimePointSpec, typename SetSpec, bool isInterval>
DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::DDEForwardTaylorCurvePiece(
		DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval> const & orig
):
	m_t0(orig.m_t0), m_dimension(orig.m_dimension), m_order(orig.m_order),
	m_jet_at_t0(0),
	m_r0(0), m_r0_owner(orig.m_r0_owner),
	m_Xi(0), m_Xi_owner(orig.m_Xi_owner)
{
	// we need to handle memory allocations carefuly
	if (orig.m_r0_owner){ allocateR0(orig.storageN0()); *m_r0 = *orig.m_r0; } else { m_r0 = orig.m_r0; }
	if (orig.m_Xi_owner){ allocateXi(); *m_Xi = *orig.m_Xi; } else { m_Xi = orig.m_Xi; }
	allocateJet();
	auto thisJetIt = beginJet();
	auto origJetIt = orig.beginJet();
	for ( ; origJetIt != orig.endJet(); ++origJetIt, ++thisJetIt) *thisJetIt = *origJetIt;
	updateCommonR0();
}

template<typename TimePointSpec, typename SetSpec, bool isInterval>
DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::DDEForwardTaylorCurvePiece(
		TimePointType t0,
		size_type dimension,
		size_type order):
	m_t0(t0), m_dimension(dimension), m_order(order),
	m_jet_at_t0(0),
	m_r0(0), m_r0_owner(true),
	m_Xi(0), m_Xi_owner(true)
{
	allocate();
	setAsConstant(VectorType(dimension));
}

template<typename TimePointSpec, typename SetSpec, bool isInterval>
DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::DDEForwardTaylorCurvePiece(
		TimePointType t0,
		size_type order,
		const VectorType& v,
		size_type N0):
	m_t0(t0), m_dimension(v.dimension()), m_order(order),
	m_jet_at_t0(0),
	m_r0(0), m_r0_owner(true),
	m_Xi(0), m_Xi_owner(true)
{
	allocateJet();
	setAsConstant(v, N0);
}

template<typename TimePointSpec, typename SetSpec, bool isInterval>
DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::DDEForwardTaylorCurvePiece(
		TimePointType t0,
		size_type order,
		SetType& v):
	m_t0(t0), m_dimension(VectorType(v).dimension()), m_order(order),
	m_jet_at_t0(0),
	m_r0(0), m_r0_owner(true),
	m_Xi(0), m_Xi_owner(true)
{
	allocateJet();
	setAsConstant(v);
}

template<typename TimePointSpec, typename SetSpec, bool isInterval>
DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::DDEForwardTaylorCurvePiece(
		TimePointType t0,
		size_type order,
		const SetType& v):
	m_t0(t0), m_dimension(VectorType(v).dimension()), m_order(order),
	m_jet_at_t0(0),
	m_r0(0), m_r0_owner(true),
	m_Xi(0), m_Xi_owner(true)
{
	allocateJet();
	setAsConstant(v);
}

template<typename TimePointSpec, typename SetSpec, bool isInterval>
DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::DDEForwardTaylorCurvePiece(
		TimePointType t0,
		size_type order,
		const VectorType& v,
		VectorType* r0,
		bool passOwnership):
	m_t0(t0), m_dimension(VectorType(v).dimension()), m_order(order),
	m_jet_at_t0(0),
	m_r0(0), m_r0_owner(true),
	m_Xi(0), m_Xi_owner(true)
{
	allocateJet();
	setAsConstant(v, r0, passOwnership);
}

template<typename TimePointSpec, typename SetSpec, bool isInterval>
DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::DDEForwardTaylorCurvePiece(
		TimePointType t0,
		VectorType const* it,
		VectorType const* itEnd):
	m_t0(t0), m_dimension(0), m_order(0),
	m_jet_at_t0(0),
	m_r0(0), m_r0_owner(true),
	m_Xi(0), m_Xi_owner(true)
{
	try { setupFromData(it, itEnd); } catch (std::logic_error& e) { throw rethrow("DDEForwardTaylorCurvePiece(t, VectorType*...): could not setup from data", e); }
}

template<typename TimePointSpec, typename SetSpec, bool isInterval>
DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::DDEForwardTaylorCurvePiece(
		TimePointType t0,
		std::vector<VectorType> coeffs):
	m_t0(t0), m_dimension(0), m_order(0),
	m_jet_at_t0(0),
	m_r0(0), m_r0_owner(true),
	m_Xi(0), m_Xi_owner(true)
{
	try { setupFromData(coeffs.begin(), coeffs.end()); } catch (std::logic_error& e) { throw rethrow("DDEForwardTaylorCurvePiece(t, std::vector... ): could not setup from data", e); }
}

template<typename TimePointSpec, typename SetSpec, bool isInterval>
DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::DDEForwardTaylorCurvePiece(
		TimePointType t0,
		std::vector<VectorType> coeffs,
		std::vector<MatrixType> Cs,
		VectorType* r0, bool passR0Ownership,
		std::vector<MatrixType> Bs,
		std::vector<VectorType> rs,
		VectorType* Xi, bool passXiOwnership):
	m_t0(t0), m_dimension(0), m_order(0),
	m_jet_at_t0(0),
	m_r0(0), m_r0_owner(true),
	m_Xi(0), m_Xi_owner(true)
{
	if (rs.size() != Bs.size())
		throw std::logic_error("DDEForwardTaylorCurvePiece::__construct__(): Bs is not equal in size to rs");
	if (rs.size() != coeffs.size())
		throw std::logic_error("DDEForwardTaylorCurvePiece::__construct__(): rs is not equal in size to coeffs");
	MatrixType Z(coeffs[0].dimension(), coeffs[0].dimension());
	std::vector<MatrixType> invBs(Bs.size(), Z);
	for (size_type k = 0; k < coeffs.size(); ++k)
		Class::Policy(Bs[k], invBs[k], rs[k]);
		//this->Policy(Bs[k], invBs[k], rs[k]);
	try {
		setupFromData(
			coeffs.begin(), coeffs.end(),
			Cs.begin(), Cs.end(),
			r0,
			Bs.begin(), Bs.end(),
			rs.begin(), rs.end(),
			Xi
		);
		m_r0_owner = passR0Ownership;
		m_Xi_owner = passXiOwnership;
	} catch (std::logic_error& e) {
		helper_safe_delete(r0, passR0Ownership);
		helper_safe_delete(Xi, passXiOwnership);
		throw rethrow("DDEForwardTaylorCurvePiece(t, full... ): could not setup from data", e);
	}
}
template<typename TimePointSpec, typename SetSpec, bool isInterval>
DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::DDEForwardTaylorCurvePiece(
		TimePointType t0,
		std::vector<VectorType> coeffs,
		std::vector<MatrixType> Cs,
		VectorType* r0, bool passR0Ownership,
		std::vector<MatrixType> Bs,
		std::vector<MatrixType> invBs,
		std::vector<VectorType> rs,
		VectorType* Xi, bool passXiOwnership):
	m_t0(t0), m_dimension(0), m_order(0),
	m_jet_at_t0(0),
	m_r0(0), m_r0_owner(true),
	m_Xi(0), m_Xi_owner(true)
{
	try {
		setupFromData(
			coeffs.begin(), coeffs.end(),
			Cs.begin(), Cs.end(),
			r0,
			Bs.begin(), Bs.end(),
			invBs.begin(), invBs.end(),
			rs.begin(), rs.end(),
			Xi
		);
		m_r0_owner = passR0Ownership;
		m_Xi_owner = passXiOwnership;
	} catch (std::logic_error& e) {
		helper_safe_delete(r0, passR0Ownership);
		helper_safe_delete(Xi, passXiOwnership);
		throw rethrow("DDEForwardTaylorCurvePiece(t, full... ): could not setup from data", e);
	}
}

template<typename TimePointSpec, typename SetSpec, bool isInterval>
DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>
DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::midCurve() const {
	typedef DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval> Class;
	typedef typename Class::DoubletonStorageType Storage;
	Class result(m_t0, dimension(), order());
	auto jetIt = beginJet();
	auto resIt = result.beginJet();
	for (; jetIt != endJet(); ++jetIt, ++resIt){
		(*resIt) = Storage(jetIt->midPoint());
	}
	*result.m_Xi = *m_Xi;
	// m_r0 stays as Vector(0) from the constructor!
	return result;
}

template<typename TimePointSpec, typename SetSpec, bool isInterval>							// template spec
DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>									// return type
DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::dt(size_type n) const {	// function decl
	if (m_order < n){
		std::ostringstream info;
		info << "DDEForwardTaylorCurvePiece::dt(n): curve order < n for n = " << n << ", curve order = " << m_order << ".";
		throw std::invalid_argument(info.str());
	}
	DDEForwardTaylorCurvePiece result(m_t0, m_dimension, m_order - n);
	// TODO: (FUTURE) implement
	throw std::logic_error("DDEForwardTaylorCurvePiece::dt: Not implemented yet");
	return result;
}

template<typename TimePointSpec, typename SetSpec, bool isInterval>									// template spec
typename DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::VectorType					// return type
DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::evalAtDelta(const RealType& delta_t) const {	// function decl
	return taylorAtDelta(delta_t) + summaAtDelta(delta_t);
}

template<typename TimePointSpec, typename SetSpec, bool isInterval>													// template spec
typename DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::VectorType									// return type
DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::evalCoeffAtDelta(size_type n, const RealType& delta_t) const {	// function decl
	if (n == 0) return this->evalAtDelta(delta_t);
	else if (n == m_order + 1) return *(this->m_Xi); 	// this is very rough, but we cannot produce better!
	else if (n > m_order + 1) {							// we cannot guarantee the continuity of higher class.
		std::ostringstream info;
		info << "DDEForwardTaylorCurvePiece::evalCoeff: cannot eval derivative of higher order. ";
		info << "Max order available: " << m_order + 1 << ", requested order: " << n << ".";
		throw std::logic_error(info.str());
	}
	// otherwise we need to do some computation...
	VectorType result = *(this->m_Xi);
	size_type j = this->m_order;
	const_iterator it = this->backJet();
	while (true){
		result = VectorType(*it) + result * (delta_t * (RealType(j + 1) / RealType(j + 1 - n)));
		if (j == n) break; else { --j; --it; } // makes sure that j is good if size_type is unsigned! (case n=0)
	}
	return result;
}

template<typename TimePointSpec, typename SetSpec, bool isInterval>
void DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::evalAtDelta(
			const RealType& delta_t, SetSpec& out) const {
	this->taylorAtDelta(delta_t, out);
	// put summa part in B*r part
	out.add(summaAtDelta(delta_t));
}

template<typename TimePointSpec, typename SetSpec, bool isInterval>
void DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::evalCoeffAtDelta(
			size_type n, const RealType& delta_t, SetSpec& out) const {
	if (out.dimension() != this->dimension())
		throw std::logic_error("DDEForwardTaylorCurvePiece::evalCoeffAtDelta(): output has different dimension;");
	if (out.storageN0() != this->storageN0())
		throw std::logic_error("DDEForwardTaylorCurvePiece::evalCoeffAtDelta(): output has different N0 dimension;");
	if (n == 0) { this->evalAtDelta(delta_t, out); return; }
	else if (n == m_order + 1) { out.add(*(this->m_Xi)); return; }
	else if (n > m_order + 1) {
		std::ostringstream info;
		info << "DDEForwardTaylorCurvePiece::evalCoeff: cannot eval derivative of higher order. ";
		info << "Max order available: " << m_order + 1 << ", requested order: " << n << ".";
		throw std::logic_error(info.str());
	}
	// Tried to do .mulThenAdd() to DRY but it produces undesirable effects than
	// applying mean value form to the evaluation (the old way).
	// mulThenAdd reorganizes the set at each step which might cause dependency error
	// (when multipled again and again by delta_t)
	size_type j = this->m_order;
	size_type d = this->dimension();
	size_type N0 = this->storageN0();
	VectorType* x = new VectorType(d);
	MatrixType* C = new MatrixType(d, N0);
	// we will forget to reorganize B * r part - we evaluate it as Interval set here.
	// it would be not feasible to maintain this structure in the algorithm.
	(*x) += *(this->m_Xi);
	const_iterator coeff = this->backJet();
	while(true) {
		(*x) = coeff->get_x() + coeff->get_B() * coeff->get_r() + (RealType(j + 1) / RealType(j + 1 - n)) * delta_t * (*x);
		(*C) = coeff->get_C() + (RealType(j + 1) / RealType(j + 1 - n)) * delta_t * (*C);
		if (j == n) break; else { --j; --coeff; }
	}
	MatrixType Id(d, d); Id.setToIdentity();
	VectorType s(d); MatrixType S(d, N0);
	capd::vectalg::split(*x, s);
	capd::vectalg::split(*C, S);
	s += S * (*m_r0);
	out.set_x(x, true);
	out.set_Cr0(C, m_r0, true, false);
	out.set_B(Id);
	out.set_r(s);
}

template<typename TimePointSpec, typename SetSpec, bool isInterval>										// template spec
typename DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::VectorType						// return type
DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::taylorAtDelta(const RealType& delta_t) const {		// function decl
	VectorType result(dimension());
	const_iterator coeff = endJet();
	const_iterator stop = beginJet();
	while(coeff-- != stop)
		result = VectorType(*coeff) + delta_t * result;

	return result;
}

template<typename TimePointSpec, typename SetSpec, bool isInterval>
void DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::taylorAtDelta(
			const RealType& delta_t, SetSpec& out) const {
	const_iterator coeff = endJet();
	const_iterator stop = beginJet();
	out *= 0.0; // this will assure ownership... and zero the output (do I want this?)
	size_type d = this->dimension();
	size_type N0 = this->storageN0();
	VectorType* x = new VectorType(d);
	MatrixType* C = new MatrixType(d, N0);
	while(coeff-- != stop){
		(*x) = coeff->get_x() + coeff->get_B() * coeff->get_r() + delta_t * (*x);
		(*C) = coeff->get_C() + delta_t * (*C);
	}

	VectorType s(d); MatrixType S(d, N0);
	MatrixType Id(d, d); Id.setToIdentity();
	// we use Id as B part, as we do not want to decide how to choose B from the taylor mess.
	capd::vectalg::split(*x, s);
	capd::vectalg::split(*C, S);
	s += S * (*m_r0);
	out.set_x(x, true);
	out.set_Cr0(C, m_r0, true, false); // TODO: (SERIOUS RETHINK): now the out and this curve piece shares r0, might be dangerous if this piece disappears (I had some issues with that)
	out.set_B(Id);
	out.set_Binv(Id);
	out.set_r(s);
}


template<typename TimePointSpec, typename SetSpec, bool isInterval>										// template spec
typename DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::VectorType						// return type
DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::summaAtDelta(const RealType& delta_t) const {		// function decl
	// NOTE: (n+1) \int_0^t \xi(s) (t-s)^n ds with \xi(s) \in M - interval, translates
	//       into bound of type M * t^{n+1}, since
	//       (n+1) \int_0^t \xi(s) (t-s)^n ds \in M * (n+1) \int_0^t (t-s)^n ds = M * (n+1)/(n+1) (t-s)^{n+1}|_0^t
	return (get_Xi() * power(delta_t, order() + 1));
}

template<typename TimePointSpec, typename SetSpec, bool isInterval>						// template spec
typename DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::VectorType		// return type
DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::get_x() const {		// function decl
	VectorType result = makeStorage_x();
	size_type I = 0;
	// TODO: (NOT URGENT) rewrite with iterators (see taylorAtDelta(const RealType& delta_t) for inspiration)
	for (size_type k = 0; k <= order(); ++k){
		VectorType x = m_jet_at_t0[k].get_x();
		for (size_type j = 0; j < dimension(); j++)
			result[I++] = x[j];
	}
	return result;
}

template<typename TimePointSpec, typename SetSpec, bool isInterval>						// template spec
typename DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::MatrixType		// return type
DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::get_C() const {		// function decl
	MatrixType result = makeStorage_C();
	size_type I = 0;
	// TODO: (NOT URGENT) rewrite with iterators (see taylorAtDelta(const RealType& delta_t) for inspiration)
	for (size_type k = 0; k <= order(); ++k){
		MatrixType C = m_jet_at_t0[k].get_C();
		for (size_type j = 0; j < dimension(); j++)
			result.row(I++) = C.row(j);
	}
	return result;
}

template<typename TimePointSpec, typename SetSpec, bool isInterval>						// template spec
typename DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::VectorType		// return type
DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::get_r0() const {		// function decl
	return *m_r0;
}

template<typename TimePointSpec, typename SetSpec, bool isInterval>						// template spec
typename DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::MatrixType		// return type
DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::get_B() const {		// function decl
	MatrixType result = makeStorage_B();
	size_type I = 0;
	// TODO: (NOT URGENT) rewrite with iterators (see taylorAtDelta(const RealType& delta_t) for inspiration)
	for (size_type k = 0; k <= order(); ++k){
		MatrixType B = m_jet_at_t0[k].get_B();
		for (size_type i = 0; i < dimension(); i++)
			for (size_type j = 0; j < dimension(); j++)
				result[I + i][I + j] = B[i][j];
		I += dimension();
	}
	return result;
}

template<typename TimePointSpec, typename SetSpec, bool isInterval>						// template spec
typename DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::VectorType		// return type
DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::get_r() const {		// function decl
	VectorType result = makeStorage_x();
	size_type I = 0;
	// TODO: (NOT URGENT) rewrite with iterators (see taylorAtDelta(const RealType& delta_t) for inspiration)
	for (size_type k = 0; k <= order(); ++k)
		for (size_type j = 0; j < dimension(); j++)
			result[I++] = m_jet_at_t0[k].get_r()[j];
	return result;
}

template<typename TimePointSpec, typename SetSpec, bool isInterval>					// template spec
typename DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::BaseClass& 	// return type
DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::set_x(				// function decl
		VectorType const &x) 																// params
{
	size_type I = 0;
	VectorType v(dimension());
	// TODO: (NOT URGENT) rewrite with iterators (see taylorAtDelta(const RealType& delta_t) for inspiration)
	for (size_type k = 0; k <= order(); ++k){
		for (size_type j = 0; j < dimension(); j++)
			v[j] = x[I++];
		m_jet_at_t0[k].set_x(v);
	}
	return *this;
}

template<typename TimePointSpec, typename SetSpec, bool isInterval>					// template spec
typename DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::BaseClass& 	// return type
DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::set_C(				// function decl
		MatrixType const &C) 																// params
{
	size_type I = 0;
	MatrixType M(dimension(), storageN0());
	// TODO: (NOT URGENT) rewrite with iterators (see taylorAtDelta(const RealType& delta_t) for inspiration)
	for (size_type k = 0; k <= order(); ++k){
		for (size_type j = 0; j < dimension(); j++)
			M.row(j) = C.row(I++);
		m_jet_at_t0[k].set_C(M);
	}
	return *this;
}

template<typename TimePointSpec, typename SetSpec, bool isInterval>					// template spec
typename DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::BaseClass& 	// return type
DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::set_Cr0(			// function decl
		MatrixType const &C, VectorType const &r0)											// params
{
	// TODO: (NOT URGENT) rewrite with iterators (see taylorAtDelta(const RealType& delta_t) for inspiration)
	if (C.numberOfColumns() != r0.dimension())
		throw std::logic_error("DDEForwardTaylorCurvePiece::set_Cr0: incompatible dimensions.");
	// this chain of instructions should work on all SetTypes
	// TODO: (???) for SharedSet it could be done easier
	*m_r0 = r0;
	size_type I = 0;
	for (size_type k = 0; k <= order(); ++k){
		MatrixType* M = new MatrixType(dimension(), storageN0());
		for (size_type j = 0; j < dimension(); j++)
			M->row(j) = C.row(I++);
		// we use set_Cr0(ptr, ptr, own, own) command in what follows
		// we pass true as the ownership of M matrix, and false (we keep ownership of m_r0)
		// it would be very fast for SharedSet (copy ptrs) but slower for BasicSet (essentially copy by value)
		m_jet_at_t0[k].set_Cr0(M, m_r0, true, false);
	}
	return *this;
}

template<typename TimePointSpec, typename SetSpec, bool isInterval>					// template spec
typename DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::BaseClass& 	// return type
DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::set_B(				// function decl
		MatrixType const &B) 																// params
{
	if (B.numberOfColumns() != B.numberOfRows())
		throw std::logic_error("DDEForwardTaylorCurvePiece::set_B(): B must be square");
	if (B.numberOfColumns() != storageDimension())
		throw std::logic_error("DDEForwardTaylorCurvePiece::set_B(): B has dimensions not compatible with the set");
	size_type offCount = 0;
	std::vector<MatrixType> Bs = extractDiagonalBlocks(B, dimension(), offCount);
	if (offCount)
		throw std::logic_error("DDEForwardTaylorCurvePiece::set_B(): off-block-diagonal element of B is nonzero. Check documentation why this is an error.");

	auto b = Bs.begin();
	auto c = this->beginJet();
	for (; b != Bs.end(); ++b, ++c)
		c->set_B(*b);
	return *this;
}

template<typename TimePointSpec, typename SetSpec, bool isInterval>					// template spec
typename DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::BaseClass& 	// return type
DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::set_Binv(			// function decl
		MatrixType const &invB) 															// params
{
	// TODO: (NOT URGENT) DRY (see set_B())
	if (invB.numberOfColumns() != invB.numberOfRows())
		throw std::logic_error("DDEForwardTaylorCurvePiece::set_Binv(): invB must be square");
	if (invB.numberOfColumns() != storageDimension())
		throw std::logic_error("DDEForwardTaylorCurvePiece::set_Binv(): invB has dimensions not compatible with the set");
	size_type offCount = 0;
	std::vector<MatrixType> Bs = extractDiagonalBlocks(invB, dimension(), offCount);
	if (offCount)
		throw std::logic_error("DDEForwardTaylorCurvePiece::set_Binv(): off-block-diagonal element of invB is nonzero. Check documentation why this is an error.");

	auto b = Bs.begin();
	auto c = this->beginJet();
	for (; b != Bs.end(); ++b, ++c)
		c->set_Binv(*b);
	return *this;
}

template<typename TimePointSpec, typename SetSpec, bool isInterval>					// template spec
typename DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::BaseClass& 	// return type
DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::set_r(				// function decl
		VectorType const &r) 																// params
{
	// TODO: (NOT URGENT) rewrite with iterators (see taylorAtDelta(const RealType& delta_t) for inspiration)
	size_type I = 0;
	VectorType v(dimension());
	for (size_type k = 0; k <= order(); ++k){
		for (size_type j = 0; j < dimension(); j++)
			v[j] = r[I++];
		m_jet_at_t0[k].set_r(v);
	}
	return *this;
}

template<typename TimePointSpec, typename SetSpec, bool isInterval>					// template spec
typename DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::BaseClass& 	// return type
DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::set_r0(			// function decl
		VectorType const &r0) 																// params
{
	// if dimension is bad, then setting r0 would destroy the sanity of the doubleton structure.
	if (r0.dimension() != storageN0())
		throw std::logic_error("DDEForwardTaylorCurvePiece::set_r0(): new r0 has incompatible N_0 dimension.");
	*m_r0 = r0; // we copy the value of the vector, so we do not need to updateR0 (pointers are still good).
	return *this;
}


template<typename TimePointSpec, typename SetSpec, bool isInterval>
std::string
DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval>::show() const {
	std::ostringstream result;
	result << "DDEForwardTaylorCurvePiece\n";
	result << "c(t) = " << m_jet_at_t0[0].hull();
	std::string prefix = "";
	if (order() > 0) result << " + " << m_jet_at_t0[1].hull() << " * (t - " << m_t0 << ")";
	for (size_type k = 2; k <= order(); ++k){
		result << " + " << m_jet_at_t0[k].hull() << " * (t - " << m_t0 << ")^" << k;
	}
	result << " + " << (order()+1) << " * \\int_{" << m_t0 << "}^{t} \\xi(s) * (t-s)^" << order() << "ds";
	result << ", with \\xi(s) \\in " << get_Xi();
	return result.str();
}

} // namespace ddes
} // namespace capd

#endif /* _CAPD_DDES_DDECURVEPIECE_HPP_ */
