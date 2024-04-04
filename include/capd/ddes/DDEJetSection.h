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

#ifndef _CAPD_DDES_DDE_JET_SECTION_H_
#define _CAPD_DDES_DDE_JET_SECTION_H_

#include <capd/ddes/DDECommon.h>
#include <capd/ddes/storage/GenericJet.h>

namespace capd{
namespace ddes{

/**
 * Represents a section that is a hyperplane with a given normal vector in some function space.
 *
 * Currently, the hyperplane is given as S = { v : <normal . v> = c }, where the scalar product < . > is
 * the standard product in \R^M, where M is the dimension of all coefficients defining the space.
 *
 * It can bu used e.g. to represent many standard
 *
 * TODO: (FUTURE) allow the Curve to define .dot() product for the hyperplane, to allow for more natural functional space products, i.e. L^2 or such.
 * TODO: (NOTE) consider e.g. L^2 product: s . v = \int_0^h s(t) * v(t) dt, where:
 *    s(t) = s0 + s1*t + s2*t^2 + ... + sn*t^n
 *    v(t) = v0 + v1*t + v2*t^2 + ... + vn*t^n
 *    s . v = \sum_j=0^n v0 * (\int_0^h h^j * s(t) dt ), the last integral depends only on s, can be computed explicitly.
 *    The j-th integral in ( ) is the j-th coefficent of the normal vector you should supply to
 *    have L^2 norm. The coefficients of course should be applied jet-wise, over intervals [-ih, -ih +h).
 */
template<typename CurveSpec, bool isInterval = capd::TypeTraits<typename CurveSpec::ScalarType>::isInterval>
class DDEJetSection {
public:
	typedef CurveSpec CurveType;
	typedef typename CurveType::size_type size_type;
	typedef typename CurveType::TimePointType TimePointType;
	typedef typename CurveType::VectorType VectorType;
	typedef typename CurveType::MatrixType MatrixType;
	typedef typename CurveType::ScalarType ScalarType;
	typedef GenericJet<TimePointType, VectorType, VectorType, MatrixType, isInterval> JetType;
	typedef std::vector<JetType> JetStorageType;
	typedef typename JetStorageType::const_iterator const_iterator;
	typedef typename JetStorageType::iterator iterator;
	typedef typename JetStorageType::const_reverse_iterator const_reverse_iterator;
	typedef typename JetStorageType::reverse_iterator reverse_iterator;

	DDEJetSection(){}

	/**
	 * makes a section checking value at t=0 of a given direction, i.e. coordinate section:
	 * x(0)_i == c, where x(0) is d dimensional.
	 */
	DDEJetSection(size_type d, size_type i, ScalarType c = 0.): m_c(c) {
		if (i < 0 || i >= d) // TODO: (NOT URGENT) more info in exception
			throw std::logic_error("DDEJetSection::__construct__(int, int, scalar): i must be between 0 and d.");
		m_orig_s = VectorType(d);
		m_orig_s[i] = 1.0;
		m_s = m_orig_s;
	}

	/**
	 * makes a section with p Jets of order n each and expected scalar value c.
	 * vector vec must be d * (1 + p * (n+1)) dimensional and will propagate the jets.
	 */
	DDEJetSection(size_type d, size_type p, size_type n, VectorType const& vec, ScalarType c): m_c(c) {
		size_type expectedDim = d * (1 + p * (n+1));
		if (expectedDim != vec.dimension()) // TODO: (NOT URGENT) more info in exception
			throw std::logic_error("DDEJetSection::__construct__(): vector has bad dimension.");

		m_orig_s = vec;
		size_type I = 0;
		m_s = VectorType(d); for (size_type i = 0; i < d; ++i, ++I) m_s[i] = vec[I];
		for (size_type j = 0; j < p; ++j)
			m_jets.push_back(JetType(TimePointType(), d, n));
		for (auto jet = rbegin(); jet != rend(); ++jet)
			for (size_type k = 0; k <= n; ++k)
				for (size_type i = 0; i < d; ++i, ++I)
					(*jet)[k][i] = vec[I];
	}

	/**
	 * makes a section that is a hyperplane with a given normal vector in some
	 * function space (where @param normal lives, i.e. C^\eta_p space).
	 *
	 * It inherits structure from the normal, and the hyperplane is
	 * S = { v : normal . v = c }, where the scalar product . is
	 * the standard product in the
	 */
	DDEJetSection(CurveSpec const& normal, ScalarType c): m_c(c) {
		VectorType vec = normal;
		size_type expectedDim = vec.dimension();

		m_orig_s = vec;
		m_s = normal.getValueAtCurrent();
		for (auto jet: normal)
			m_jets.push_back(*jet);
	}

	ScalarType operator()(CurveSpec const& curve) const {
		// TODO: (NOT URGENT) dim checking?
		return curve.dot(*this) - m_c;
	}

	/**
	 * This returns gradient of the section d/dx s(x) evaluated at x = j(curve)
	 *
	 * NOTE: Since in our case section is a hypersurface, the gradient is constant
	 *       and equal (z(s), j(s)), extended to the structure of curve, se next NOTE.
	 *
	 * NOTE: this function uses structure in curve to deliver higher order
	 *       jet elements needed for computations. They will be set simply to 0.
	 *       Therefore, it would be possible to multiply as a Vectors, i.e.
	 *       we can simply do the vector-vector dot product:
	 *
	 *       	getGradient(curve) . (z(curve), j(curve))
	 *
	 *       without worrying about (z(s), j(s)) having bad dimensions.
	 *
	 * Note: curve passed by copy, as I need it inside
	 */
	VectorType getGradient(CurveSpec curve) const {
		curve *= 0.; // i need only the structure!
		auto icjet = curve.rbegin();
		auto isjet = rbegin();
		for (; icjet != curve.rend() && isjet != rend(); ++icjet, ++isjet){
			auto iccoeff = (*icjet)->begin();
			auto iscoeff = isjet->begin();
			for (; iscoeff != isjet->end(); ++iscoeff, ++iccoeff)
				(*iccoeff) = (*iscoeff);
		}
		curve.setValueAtCurrent(m_s);
		return curve.get_x();
	}

	operator VectorType() const {
		return this->m_orig_s;
	}

	size_type storageDimension() const { return m_orig_s.dimension(); }
	size_type dimension() const { return m_s.dimension(); }

	// TODO: (NOT URGENT) setters, getters, extenders.

	DDEJetSection& set_c(ScalarType c){ m_c = c; return *this; }
	ScalarType get_c() const { return m_c; }
	DDEJetSection& set_s(VectorType s){
		if (s.dimension() != m_s.dimension()) // TODO: (NOT URGENT) more info in exception
			throw std::logic_error("DDEJetSection::set_s(): vector has bad dimension.");
		// TODO: (NOT URGENT) check if not contains 0. vector (?), TODO: (NOT URGENT) check dimension?
		for (size_type i = 0; i < m_s.dimension(); ++i) m_orig_s[i] = s[i];
		m_s = s;
		return *this;
	}
	VectorType const& get_s() const { return m_s; }
	DDEJetSection& extend(JetType const& jet){
		if (jet.dimension() != m_s.dimension()) // TODO: (NOT URGENT) more info in exception
			throw std::logic_error("DDEJetSection::extend(): jet has bad dimension.");
		m_jets.push_back(jet);
		size_type n = jet.order();
		size_type d = m_s.dimension();
		VectorType new_orig_s(m_orig_s.dimension() + d * (n+1));
		size_type i = 0; for (; i < m_orig_s.dimension(); ++i) new_orig_s[i] = m_orig_s[i];
		for (auto coeff = jet.begin(); coeff != jet.end(); ++coeff)
			for (size_type j = 0; j < d; ++j)
				new_orig_s[i++] = (*coeff)[j];
		m_orig_s = new_orig_s;
		return *this;
	}

	reverse_iterator rbegin(){ return m_jets.rbegin(); }
	const_reverse_iterator rbegin() const { return m_jets.rbegin(); }
	reverse_iterator rend(){ return m_jets.rend(); }
	const_reverse_iterator rend() const { return m_jets.rend(); }
	iterator begin(){ return m_jets.begin(); }
	const_iterator begin() const { return m_jets.begin(); }
	iterator end(){ return m_jets.end(); }
	const_iterator end() const { return m_jets.end(); }


protected:
	ScalarType m_c;
	VectorType m_s;
	JetStorageType m_jets;

	VectorType m_orig_s;
};

} // namespace ddes
} // namespace capd

#endif /* _CAPD_DDES_DDE_JET_SECTION_H_ */
