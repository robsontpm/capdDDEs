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

#ifndef _CAPD_DDES_GENERICJET_HPP_
#define _CAPD_DDES_GENERICJET_HPP_

#include <capd/ddes/DDECommon.hpp>
#include <capd/ddes/storage/GenericJet.h>

////////////////////////////////////////////////////////////////////////
// Below are implementations of more involved functions               //
// (to be separated as usual from the declarations in the .h file     //
// For more complete documentation/explanation please consult .h file //
////////////////////////////////////////////////////////////////////////

namespace capd{
namespace ddes{

template<typename TimePointSpec, typename DataSpec, typename VectorSpec, typename MatrixSpec, bool isInterval>
GenericJet<TimePointSpec, DataSpec, VectorSpec, MatrixSpec, isInterval>&
GenericJet<TimePointSpec, DataSpec, VectorSpec, MatrixSpec, isInterval>::operator=(
		GenericJet<TimePointSpec, DataSpec, VectorSpec, MatrixSpec, isInterval> const & orig
)
{
	m_dimension = orig.m_dimension;
	m_order = orig.m_order;
	m_t0 = orig.m_t0;
	setupCoeffs(orig.begin(), orig.end()); // not need to try catch, orig assumed to be sane.
	return *this;
}

template<typename TimePointSpec, typename DataSpec, typename VectorSpec, typename MatrixSpec, bool isInterval>
GenericJet<TimePointSpec, DataSpec, VectorSpec, MatrixSpec, isInterval>::GenericJet(
		GenericJet<TimePointSpec, DataSpec, VectorSpec, MatrixSpec, isInterval> const & orig
):
	m_t0(orig.m_t0), m_dimension(orig.m_dimension), m_order(orig.m_order), m_coeffs(0)
{
	setupCoeffs(orig.begin(), orig.end()); // not need to try catch, orig assumed to be sane.
}

template<typename TimePointSpec, typename DataSpec, typename VectorSpec, typename MatrixSpec, bool isInterval>
GenericJet<TimePointSpec, DataSpec, VectorSpec, MatrixSpec, isInterval>::GenericJet(
		TimePointType t0,
		size_type dimension,
		size_type order):
	m_t0(t0), m_dimension(dimension), m_order(order), m_coeffs(0)
{
	allocateCoeffs();
}

template<typename TimePointSpec, typename DataSpec, typename VectorSpec, typename MatrixSpec, bool isInterval>
GenericJet<TimePointSpec, DataSpec, VectorSpec, MatrixSpec, isInterval>::GenericJet(
		TimePointType t0,
		size_type order,
		const DataType& v):
	m_t0(t0), m_dimension(v.dimension()), m_order(order), m_coeffs(0)
{
	allocateCoeffs();
	m_coeffs[0] = v;
}

template<typename TimePointSpec, typename DataSpec, typename VectorSpec, typename MatrixSpec, bool isInterval>
GenericJet<TimePointSpec, DataSpec, VectorSpec, MatrixSpec, isInterval>::GenericJet(
		TimePointType t0,
		const_iterator it,
		const_iterator itEnd):
	m_t0(t0), m_dimension(0), m_order(0), m_coeffs(0)
{
	try { setupCoeffs(it, itEnd); }
	catch (std::logic_error& e) { deallocateCoeffs(); throw rethrow("GenericJet(t, VectorType*...): could not setup from data", e); }
}

template<typename TimePointSpec, typename DataSpec, typename VectorSpec, typename MatrixSpec, bool isInterval>
GenericJet<TimePointSpec, DataSpec, VectorSpec, MatrixSpec, isInterval>::GenericJet(
		TimePointType t0,
		size_type order,
		const_iterator it,
		const_iterator itEnd):
	m_t0(t0), m_dimension(0), m_order(0), m_coeffs(0)
{
	try { setupCoeffs(it, itEnd, order); }
	catch (std::logic_error& e) { deallocateCoeffs(); throw rethrow("GenericJet(t, order, VectorType*...): could not setup from data", e); }
}

template<typename TimePointSpec, typename DataSpec, typename VectorSpec, typename MatrixSpec, bool isInterval>
GenericJet<TimePointSpec, DataSpec, VectorSpec, MatrixSpec, isInterval>::GenericJet(
		TimePointType t0,
		std::vector<DataType> coeffs):
	m_t0(t0), m_dimension(0), m_order(0), m_coeffs(0)
{
	try { setupCoeffs(coeffs); }
	catch (std::logic_error& e) { deallocateCoeffs(); throw rethrow("GenericJet(t, std::vector... ): could not setup from data", e); }
}

template<typename TimePointSpec, typename DataSpec, typename VectorSpec, typename MatrixSpec, bool isInterval>	// template spec
GenericJet<TimePointSpec, DataSpec, VectorSpec, MatrixSpec, isInterval>&										// return type
GenericJet<TimePointSpec, DataSpec, VectorSpec, MatrixSpec, isInterval>::setupCoeffs(							// func decl
		const_iterator it, const_iterator itEnd, size_type overrideOrder										// params
){
	if (it == NULL || itEnd == NULL || itEnd <= it)
		throw std::logic_error("GenericJet::setupCoeffs(begin*, end*): no data given (either one of ptr is NULL or begin >= end)");
	deallocateCoeffs();
	m_dimension = it->dimension();
	m_order = overrideOrder ? overrideOrder : (itEnd - it - 1);
	allocateCoeffs();
	auto ic = this->begin();
	for (size_type i = 0; it != itEnd && i <= m_order; ++it, ++ic, ++i){
		if (it->dimension() != m_dimension){
			std::ostringstream info;
			info << "GenericJet::setupCoeffs(begin*, end*): incompatible dimension of coeff, ";
			info << "is " << it->dimension() << ", should be " << m_dimension << ", ";
			info << "at element index " << i << ".";
			deallocateCoeffs();
			throw std::logic_error(info.str());
		}
		*ic = *it;
	}
	return *this;
}

template<typename TimePointSpec, typename DataSpec, typename VectorSpec, typename MatrixSpec, bool isInterval>	// template spec
GenericJet<TimePointSpec, DataSpec, VectorSpec, MatrixSpec, isInterval>&												// return type
GenericJet<TimePointSpec, DataSpec, VectorSpec, MatrixSpec, isInterval>::setupCoeffs(							// func decl
		std::vector<DataType> const& coeffs, size_type overrideOrder														// params
){
	if (coeffs.size() == 0) throw std::logic_error("GenericJet::setupCoeffs(std::vector): no data given");
	deallocateCoeffs();
	m_dimension = coeffs[0].dimension();
	m_order = overrideOrder ? overrideOrder : (coeffs.size() - 1);
	allocateCoeffs();
	auto ic = this->begin();
	auto it = coeffs.begin();
	for (size_type i = 0; it != coeffs.end() && i <= m_order; ++it, ++ic, ++i){
		if (it->dimension() != m_dimension){
			std::ostringstream info;
			info << "GenericJet::setupCoeffs(std::vector): incompatible dimension of coeff, ";
			info << "is " << it->dimension() << ", should be " << m_dimension << ", ";
			info << "at element index " << i << ".";
			deallocateCoeffs();
			throw std::logic_error(info.str());
		}
		*ic = *it;
	}
	return *this;
}

template<typename TimePointSpec, typename DataSpec, typename VectorSpec, typename MatrixSpec, bool isInterval>	// template spec
VectorSpec																													// return type
GenericJet<TimePointSpec, DataSpec, VectorSpec, MatrixSpec, isInterval>::evalAtDelta(							// func decl
		const RealType& delta_t																								// params
) const {
	VectorType result(dimension());
	auto coeff = this->end();
	while(coeff-- != this->begin()) {
		VectorType w = *coeff;
		result = w + delta_t * result;
	}
	return result;
}

template<typename TimePointSpec, typename DataSpec, typename VectorSpec, typename MatrixSpec, bool isInterval>	// template spec
void																														// return type
GenericJet<TimePointSpec, DataSpec, VectorSpec, MatrixSpec, isInterval>::evalAtDelta(							// func decl
			const RealType& delta_t, DataSpec& out) const {
	// naive implementation, but safe for all types
	// out = DataSpec(this->evalAtDelta(delta_t));
	const_iterator coeff = this->end();
	while (coeff-- != begin()){
		out *= delta_t;   // this way it is easier to impelement neccessary (unary) operators, and no copy operator is involved! TODO: (make the same in the rigorous code)
		out += (*coeff);  // this way it is easier to impelement neccessary (unary) operators, and no copy operator is involved! TODO: (make the same in the rigorous code)
		// out = (*coeff) + out * delta_t; // old version, one line.
	}
}

template<typename TimePointSpec, typename DataSpec, typename VectorSpec, typename MatrixSpec, bool isInterval>	// template spec
VectorSpec																													// return type
GenericJet<TimePointSpec, DataSpec, VectorSpec, MatrixSpec, isInterval>::evalCoeffAtDelta(						// func decl
		size_type n, const RealType& delta_t																				// params
) const {
	if (n == 0) return this->evalAtDelta(delta_t);
	else if (n > m_order) {
		return VectorType(dimension()); // is zero
	}
	// otherwise we need to do some computation...
	VectorType result(dimension());
	size_type j = this->m_order;
	const_iterator coeff = this->end();
	while (coeff-- != begin()){
		// TODO: (tabularize the weights for some minor speed?)
		result = VectorType(*coeff) + result * (delta_t * (RealType(j + 1) / RealType(j + 1 - n)));
		if (j == n) break; // this will always be valid even if size_type is unsigned and even if n = 0
		--j;
	}
	return result;
}

template<typename TimePointSpec, typename DataSpec, typename VectorSpec, typename MatrixSpec, bool isInterval>	// template spec
void																														// return type
GenericJet<TimePointSpec, DataSpec, VectorSpec, MatrixSpec, isInterval>::evalCoeffAtDelta(						// func decl
			size_type n, const RealType& delta_t, DataSpec& out) const {
	// naive implementation, but safe for all types
	// out = DataSpec(this->evalCoeffAtDelta(n, delta_t));

	// NOTE: we assume out is of good shape (i.e. dimension, other data if neccessary), and set to 0
	if (n == 0) this->evalAtDelta(delta_t, out);
	else if (n > m_order) {
		return;
	}
	size_type j = this->m_order;
	const_iterator coeff = this->end();
	while (coeff-- != begin()){
		// TODO: (NOT URGENT, FUTURE, RETHINK) tabularize the weights for some minor speed?
		// old version in one line: out = (*coeff) + out * (delta_t * (RealType(j + 1) / RealType(j + 1 - n)));
		// is not so good as the new version, that does not use cause copy constructors, etc.
		// (just with the unary +=, *= operators):
		out *= (delta_t * (RealType(j + 1) / RealType(j + 1 - n)));
		out += (*coeff);
		// to prevent infinite loop in case of unsigned data.
		// This will always be valid even if size_type is unsigned and even if n = 0
		if (j == n) break; --j;
	}
}

template<typename TimePointSpec, typename DataSpec, typename VectorSpec, typename MatrixSpec, bool isInterval>	// template spec
std::string																													// return type
GenericJet<TimePointSpec, DataSpec, VectorSpec, MatrixSpec, isInterval>::show() const {
	std::ostringstream result;
	result << this->badge() << "\n";
	result << "c(t) = " << VectorType(m_coeffs[0]);
	std::string prefix = "";
	if (order() > 0) result << " + " << VectorType(m_coeffs[1]) << " * (t - " << RealType(m_t0) << ")";
	for (size_type k = 2; k <= order(); ++k){
		result << " + " << VectorType(m_coeffs[k]) << " * (t - " << RealType(m_t0) << ")^" << k;
	}
	return result.str();
}


} // namespace ddes
} // namespace capd

#endif /* _CAPD_DDES_GENERICJET_HPP_ */
