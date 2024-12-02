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

#ifndef _CAPD_DDES_GENERICJET_H_
#define _CAPD_DDES_GENERICJET_H_

#include <capd/ddes/DDECommon.h>
#include <capd/basicalg/TypeTraits.h>
#include <stdexcept>
#include <string>

namespace capd{
namespace ddes{

/**
 * A class to represent a generic jet for any data type
 *
 * We assume that DataType has default constructor and
 * is convertible form/to VectorType, i.e. DataType has:
 *  - a constructor DataType(VectorType) or DataType(const VectorType&)
 *  - a cast operator VectorType()
 *  - a copy operator and copy constructor.
 *  - default constructor
 *
 * Examples: capd::Vector and Sets types are compatible. ddes/storage/doubleton types are compatible.
 */
template<
	typename TimePointSpec,
	typename DataSpec,
	typename VectorSpec = typename DataSpec::VectorType,
	typename MatrixSpec = typename DataSpec::MatrixType,
	bool isInterval = capd::TypeTraits<typename MatrixSpec::ScalarType>::isInterval
>
class GenericJet {
public:
	typedef TimePointSpec TimePointType;
	typedef DataSpec DataType;
	typedef VectorSpec VectorType;
	typedef MatrixSpec MatrixType;
	typedef GenericJet<TimePointType, DataType, VectorType, MatrixType> Class;
	typedef Class BaseClass;
	typedef typename MatrixType::ScalarType ScalarType;
	typedef typename MatrixType::size_type size_type;
	typedef ScalarType RealType;  // TODO: (NOT URGENT, FAR FUTURE) eventually do TypeTraits like Daniel's code
	typedef Class JetType;

    template<
		typename OtherTimePointSpec,
		typename OtherDataSpec,
		typename OtherVectorSpec = VectorType,
		typename OtherMatrixSpec = MatrixType,
		bool OtherIsInterval = capd::TypeTraits<typename OtherMatrixSpec::ScalarType>::isInterval>
	struct rebind { typedef GenericJet<OtherTimePointSpec, OtherDataSpec, OtherVectorSpec, OtherMatrixSpec, OtherIsInterval> other; };

	typedef DataType const* const_iterator;
	typedef DataType* iterator;

	virtual GenericJet& setupCoeffs(const_iterator it, const_iterator itEnd, size_type overrideOrder = 0);
	virtual GenericJet& setupCoeffs(std::vector<DataType> const& coeffs, size_type overrideOrder = 0);

	/** copy operator */
	GenericJet& operator=(GenericJet const & orig);
	/** copy constructor */
	GenericJet(GenericJet const & orig);
	/** makes a constant function of the vector value 0 */
	GenericJet(TimePointType t0 = TimePointType(), size_type dimension = 0, size_type order = 0);
	/** makes a constant function of the vector value v */
	GenericJet(TimePointType t0, size_type order, const DataType& v);
	/** fills vector with data. Determine order from the size of data. */
	GenericJet(TimePointType t0, const_iterator it, const_iterator itEnd);
	/** fills vector with data, but sets the order to n (fills higher orders with 0 vectors) */
	GenericJet(TimePointType t0, size_type order, const_iterator it, const_iterator itEnd);
	/** fills vector with data. Determine order from the size of data. */
	GenericJet(TimePointType t0, std::vector<DataType> coeffs);

	/** dimension of the value space of the Jet (i.e. $J : \R \to \R^{dimension}$) */
	size_type dimension() const { return m_dimension; }
	/** order of the Taylor part of the Jet */
	size_type order() const { return m_order; }
	/** how many data need to store everything */
	size_type storageDimension() const { return (m_order + 1) * m_dimension; }

	/** evaluates the the jet at t0 + delta_t */
	virtual VectorType evalAtDelta(const RealType& delta_t) const;
	/** evaluates the jet at t0 + delta_t */
	virtual void evalAtDelta(const RealType& delta_t, DataType& out) const;
	/** evaluates the the jet at t */
	virtual VectorType eval(const RealType& t) const { return this->evalAtDelta(t - RealType(t0())); }
	/** evaluates the jet at t */
	virtual void eval(const RealType& t, DataType& out) const { this->evalAtDelta(t - RealType(t0()), out); }

	/** evaluates the n-th Coeff of jet at t0 + delta_t ( n-th coeff is 1/n! * jet^{(n)}(t) = jet^{[n]}(t) ) */
	virtual VectorType evalCoeffAtDelta(size_type n, const RealType& delta_t) const;
	/** evaluates the n-th time derivative of jet at t0 + delta_t ( n-th coeff is 1/n! * jet^{(n)}(t) = jet^{[n]}(t) ) */
	virtual void evalCoeffAtDelta(size_type n, const RealType& delta_t, DataType& out) const;
	/** evaluates the n-th time derivative of the jet at t ( n-th coeff is 1/n! * jet^{(n)}(t) = jet^{[n]}(t) ) */
	virtual VectorType evalCoeff(size_type n, const RealType& t) const { return this->evalCoeffAtDelta(n, t - RealType(t0())); }
	/** evaluates the n-th time derivative of jet at t ( n-th coeff is 1/n! * jet^{(n)}(t) = jet^{[n]}(t) ) */
	virtual void evalCoeff(size_type n, const RealType& t, DataType& out) const { this->evalCoeffAtDelta(n, t - RealType(t0()), out); }
	/** makes a Jet that represents n-th derivative of the jet, by default n=1 */
	GenericJet dt(size_type n = 1){
		try {
			GenericJet deriv(t0(), begin() + n, end());
			// now, rescale coeffs respectively
			for (size_type k = 1; k <= n; ++k){
				ScalarType p = k;
				for (auto icoeff = deriv.begin(); icoeff != deriv.end(); ++icoeff, p += 1.0)
					(*icoeff) *= p;
			}
			return deriv;
		} catch (std::logic_error& e) {
			std::ostringstream msg;
			msg << "GenericJet::dt(" << n << "):";
			msg << "at t = " << this->getT0() << ": could not create derivative representation!";
			throw rethrow(msg.str(), e);
		}
	}
	/** makes a Jet that represents n-th coefficient (i.e. jet^{(n)}/n! of the jet as a function of t, by default n=1 */
	GenericJet coeff(size_type n = 1){
		// TODO: SOMEHOW IMPORTANT: Not Implemented Yet
		throw std::logic_error("GenericJet::coeff(): Not Implemented Yet");
//		GenericJet deriv(t0(), begin() + n, end());
//		// now, rescale coeffs respectively
//		for (size_type k = 1; k <= n; ++k){
//			ScalarType p = k;
//			for (auto icoeff = deriv.begin(); icoeff != deriv.end(); ++icoeff, p += 1.0)
//				(*icoeff) *= p;
//		}
//		return deriv;
	}

	/** more verbose output, more human-friendly. Not rigorous-friendly. */
	virtual std::string show() const;

	/**
	 * return VectorSpec representation of the jet. It is vector of dimension storageDimension
	 * and d-dim blocks of the natural order (0-th coefficient first).
	 */
	operator VectorType() const {
		VectorType result(storageDimension());
		auto ires = result.begin();
		for (auto icoef = begin(); icoef < end(); ++icoef){
			VectorType v = VectorType(*icoef); // assume DataType is convertible to vector
			for (auto iv = v.begin(); iv != v.end(); ++iv, ++ires)
				*ires = *iv;
		}
		return result;
	}
	/** return vector representation of the jet. Other name (more explicit for VectorType), can be overwritten (virtual). */
	virtual VectorType hull() const { return *this; }

	/** set the coeficients from a vector. get/set_x() are compatible. */
	Class& set_x(VectorType const& x) {
		auto r = x.begin();
		VectorType w(dimension());
		for (auto coeff = begin(); coeff != end(); ++coeff){
			for (size_type i = 0; i < dimension(); ++i, ++r)
				w[i] = *r;
			(*coeff) = DataType(w);
		}
		return *this;
	}
	/** gets the coefficients as a single vector of dimension storageDimension() */
	VectorType get_x() const {
		VectorType result(storageDimension());
		auto r = result.begin();
		for (auto coeff = begin(); coeff != end(); ++coeff){
			VectorType w(*coeff);
			for (size_type i = 0; i < dimension(); ++i, ++r)
				*r = w[i];
		}
		return result;
	}

	/** iterator to the 0-th order Taylor coefficient (might by many-dimensions) */
	iterator begin() { return m_coeffs; }
	/** iterator to the past the order place. Standard thing in C++ */
	iterator end() { return m_coeffs + m_order + 1; }
	/** iterator to the last element in coefficients. Standard thing in C++ */
	iterator back() { return m_coeffs + m_order; }
	/** iterator to the k-th element in coefficients. Standard thing in C++ */
	iterator at(size_type k) { return m_coeffs + k; }
	/** iterator to the 0-th order Taylor coefficient (might by many-dimensions) */
	const_iterator begin() const { return m_coeffs; }
	/** iterator to the past the order place. Standard thing in C++ */
	const_iterator end() const { return m_coeffs + m_order + 1; }
	/** iterator to the last element in coefficients. Standard thing in C++ */
	const_iterator back() const { return m_coeffs + m_order; }
	/** iterator to the k-th element in coefficients. Standard thing in C++ */
	const_iterator at(size_type k) const { return m_coeffs + k; }
	/** returns n-th order Taylor coefficient. */
	DataType& operator[](size_type k){ return m_coeffs[k]; }
	/** returns n-th order Taylor coefficient. */
	DataType const& operator[](size_type k) const { return m_coeffs[k]; }

	/** set the base time of the jet */
	Class& setT0(TimePointType const & t0){ m_t0 = t0; return *this; }
	/** returns time at which this jet is located */
	TimePointType getT0() const { return m_t0; }
	/** naming convention. See getT0() */
	TimePointType t0() const { return m_t0; }

	/** TODO: docs */
	DataType* getCoeffs(){ return m_coeffs; }
	/** TODO: docs */
	Class setCoeffs(DataType* coeffs, size_type order){ deallocateCoeffs(); m_order = order; m_coeffs = coeffs; return *this; }
	/** Badge must be a single word! */
	static std::string badge() { return "GenericJet"; }

	/** TODO: docs */
	friend std::ostream& operator<<(std::ostream & out, GenericJet const & jet) {
		out << GenericJet::badge() << std::endl;
		out << jet.m_dimension << " " << jet.m_order << std::endl;
		out << jet.getT0() << std::endl;
		for (auto kth = jet.begin(); kth != jet.end(); ++kth)
			out << (*kth) << std::endl;
		return out;
	}
	/** TODO: docs */
	friend std::istream& operator>>(std::istream & in, GenericJet & jet) {
		helper_dump_badge(in);
		jet.deallocateCoeffs();
		in >> jet.m_dimension >> jet.m_order;
		jet.allocateCoeffs();
		in >> jet.m_t0;
		for (size_type k = 0; k <= jet.m_order; ++k){
			in >> jet.m_coeffs[k];
		}
		return in;
	}
	/**
	 * It return true only if the time is the same, it comes from the same grid and
	 * the order is the  same and all the coefficients match.
	 */
	bool operator==(Class const& that) {
		if (m_t0 != that.m_t0) return false;
		if (m_order != that.m_order) return false;
		auto this_kth = begin();
		auto that_kth = that.begin();
		for (; this_kth != end(); ++this_kth, ++that_kth)
			if (*this_kth != *that_kth)
				return false;
		return true;
	}
	/** see operator== */
	bool operator!=(Class const& that) {
		return !this->operator==(that);
	}
	/** TODO: docs */
	Class setOrder(size_type n){
		return Class(t0(), n, begin(), begin() + 1 + (n > m_order ? m_order : n));
	}
	/** TODO: docs */
	Class increasedOrder(size_type n = 1){
		return n == 0 ? *this : setOrder(order() + n);
	}
	/** TODO: docs */
	Class decreasedOrder(size_type n = 1){
		if (n >= m_order) throw std::logic_error("GenericJet::decreasedOrder(): cannot decrease beyond order 0.");
		return n == 0 ? *this : setOrder(order() - n);
	}

	/** standard thing */
	virtual ~GenericJet() { this->deallocateCoeffs(); };

protected:
	TimePointType m_t0;
	size_type m_dimension;
	size_type m_order;
	DataType* m_coeffs;

	/** TODO: docs */
	void allocateCoeffs(){
		m_coeffs = new DataType[m_order+1];
		for (auto ic = begin(); ic != end(); ++ic)
			*ic = DataType(VectorType(dimension()));
		// we do above construction to be compatible with eventual DataType that do not has constructor DataType(int)
		// we assume that DataType is convertible form/to VectorType
	}
	/** TODO: docs */
	void deallocateCoeffs(){ capd::ddes::helper_safe_array_delete(m_coeffs, true); }
	/** TODO: docs */
	void reallocateCoeffs(){ deallocateCoeffs(); allocateCoeffs(); }

};


} // namespace ddes
} // namespace capd


#endif /* _CAPD_DDES_GENERICJET_H_ */
