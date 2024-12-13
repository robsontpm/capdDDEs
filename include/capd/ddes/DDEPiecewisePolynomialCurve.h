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

#ifndef _CAPD_DDES_DDEPIECEWISEPOLYNOMIALCURVE_H_
#define _CAPD_DDES_DDEPIECEWISEPOLYNOMIALCURVE_H_

#include <capd/basicalg/TypeTraits.h>
#include <stdexcept>
#include <string>
#include <deque>
#include <capd/ddes/DDECommon.h>
#include <capd/ddes/storage/GenericJet.h>

namespace capd{
namespace ddes{

/**
 * This class represents a piecewise polynomial function of time
 * $x : \R \to \R^d$. The function must be over some uniform grid
 * $t_i = ih$, and the function is a polynomial of some order $n_i$
 * on $I_i = [t_i, t_{i+1})$.
 *
 * The best way to get this class is to use @see capd::ddes::NonrigorousSetup:
 *
 * typedef capd::ddes::NonrigorousSetup<MyEquation> DSetup;
 * DSetup::Solution
 *
 *
 * if you have a DDEPiecewisePolynomialCurve variable constructed like that:
 *
 * Vector v =
 * DDEPiecewisePolynomialCurve x(grid, grid(i), grid(j), v); // check available
 *
 */
template<typename GridSpec, typename JetSpec>
class DDEPiecewisePolynomialCurve {
public:
	typedef GridSpec GridType;
	typedef JetSpec JetType;
	typedef JetType CurvePieceType;
	typedef DDEPiecewisePolynomialCurve<GridType, JetType> Class;
	typedef typename JetType::DataType DataType;
	typedef typename JetType::MatrixType MatrixType;
	typedef typename JetType::VectorType VectorType;
	typedef typename JetType::ScalarType ScalarType;
	typedef ScalarType RealType; // TODO: (FAR FUTURE, NOT URGENT) do as in capd by traits
	typedef typename GridType::TimePointType TimePointType;
	typedef typename JetType::size_type size_type;
	typedef std::deque<CurvePieceType*> PiecesStorageType;
	typedef typename PiecesStorageType::iterator iterator;
	typedef typename PiecesStorageType::const_iterator const_iterator;
	typedef typename PiecesStorageType::reverse_iterator reverse_iterator;
	typedef typename PiecesStorageType::const_reverse_iterator const_reverse_iterator;

    template<typename OtherJetTypeSpec>
	struct rebind { typedef DDEPiecewisePolynomialCurve<GridSpec, OtherJetTypeSpec> other; };

	/** for identifying in exceptions and output eventual derived classes */
	virtual std::string badge() const { return "DDEPiecewisePolynomialCurve"; }

	/** copy constructor, standard thing */
	DDEPiecewisePolynomialCurve(const Class& orig):
			m_grid(orig.m_grid), m_t_current(orig.m_t_current),
			m_valueAtCurrent(orig.m_valueAtCurrent),
			m_dimension(orig.m_dimension)
	{
		copyPieces(orig.begin(), orig.end());
	}

	/** assign operator, standard thing */
	Class& operator=(const Class& orig){
		if (&m_grid != &(orig.m_grid))
			throw std::logic_error(badge() + "::operator=(orig): cannot assign solution over a different grid.");
		clear();
		m_t_current = orig.m_t_current;
		m_dimension = orig.m_dimension;
		m_valueAtCurrent = orig.m_valueAtCurrent;
		copyPieces(orig.begin(), orig.end());
		return *this;
	}
	/**
	 * It return true only if the domain is the same,
	 * the grid is the same and all the jets are the same
	 * in the sense of operator== for Jets, see the docs for them.
	 * (but basically it means the jets at the same grid points must
	 * be of the same order)
	 */
	bool operator==(Class const& that) {
		if (m_valueAtCurrent != that.m_valueAtCurrent) return false;
		if (m_t_current != that.m_t_current) return false;
		if (m_grid != that.m_grid) return false;
		if (m_pieces.size() != that.m_pieces.size()) return false;
		auto this_jet = begin();
		auto that_jet = that.begin();
		for (; this_jet != end(); ++this_jet, ++that_jet)
			if (**this_jet != (**that_jet))
				return false;
		return true;
	}
	/** see operator== */
	bool operator!=(Class const& that) {
		return !this->operator==(that);
	}

	/** makes a Curve that is d dimensional, over given grid, located initially at TimePoint 0 */
	DDEPiecewisePolynomialCurve(const GridType& grid, size_type d = 0):
			m_grid(grid), m_t_current(grid.point(0)),
			m_valueAtCurrent(VectorType(d)),
			m_dimension(d)
	{ }

	/** makes a Curve that is d dimensional, over given grid starting at a given TimePoint */
	DDEPiecewisePolynomialCurve(const GridType& grid, const TimePointType& t0, size_type d = 0):
			m_grid(grid), m_t_current(t0),
			m_valueAtCurrent(VectorType(d)),
			m_dimension(d)
	{ }

	/** creates a Curve that is a single point of given value at time point t0 */
	DDEPiecewisePolynomialCurve(const GridType& grid, const TimePointType& t0, const DataType& value):
			m_grid(grid), m_t_current(t0),
			m_valueAtCurrent(value),
			m_dimension(value.dimension())
	{ }

	/**
	 * creates a constant function with a given value on interval [t1, t0]
	 * and with given order (useful for future use as a a set in integration).
	 */
	DDEPiecewisePolynomialCurve(
				const GridType& grid,
				const TimePointType& t0,
				const TimePointType& t1,
				size_type order,
				const DataType& value):
			m_grid(grid), m_t_current(t0),
			m_valueAtCurrent(value),
			m_dimension(value.dimension())
	{
		for (TimePointType t = t0; t < t1; ++t)
			addPiece(new CurvePieceType(t, order, m_valueAtCurrent));
	}

	/**
	 * creates a constant function with a given value on interval [-p * grid(1), grid(0)]
	 * and with given order (useful for future use as a a set in integration).
	 * It is especially useful if grid was made with h = tau/p, to setup initial
	 * segment on [-tau, 0]
	 */
	DDEPiecewisePolynomialCurve(
				const GridType& grid,
				const int p,
				size_type order,
				const DataType& value):
			m_grid(grid), m_t_current(grid(-p)),
			m_valueAtCurrent(value),
			m_dimension(value.dimension())
	{
		for (TimePointType t = grid(-p); t < grid(0); ++t)
			addPiece(new CurvePieceType(t, order, m_valueAtCurrent));
	}

	/**
	 * creates a representation of the function f over the interval [t0,t1]
	 *
	 * TODO: move impl to .hpp
	 */
	template<typename AnyMatrixSpec>
	DDEPiecewisePolynomialCurve(
				const GridType& grid,
				const TimePointType& t0,
				const TimePointType& t1,
				size_type order,
				const capd::map::Map<AnyMatrixSpec>& f):
			m_grid(grid), m_t_current(t0),
			m_valueAtCurrent(VectorType(f.imageDimension())),
			m_dimension(f.imageDimension())
	{
		typedef capd::map::Map<AnyMatrixSpec> FMap;
		typedef typename FMap::JetType FJet;
		if (f.dimension() != 1){
			std::ostringstream info;
			info << "DDESolutionCurve::__construct__(): f should be R -> R^n, ";
			info << "is R^" << f.dimension() << " -> R^" << f.imageDimension() << ".";
			throw std::logic_error(info.str());
		}
		try {
			for (TimePointType t = t0; t < t1; t += grid.point(1)){
				FJet tt(f.imageDimension(), 1, order + 1);
				auto tti = tt.begin(0);
				*(tti) = ScalarType(t);
				*(++tti) = 1.;
				FJet ft = f(tt);
				std::vector<VectorType> items(order + 1, VectorType(m_dimension));
				for (size_type i = 0; i < f.imageDimension(); ++i){
					size_type j = 0;
					for (auto fti = ft.begin(i); j <= order; ++fti, ++j){
						items[j][i] = (*fti);
					}
				}
				CurvePieceType* piece = new CurvePieceType(t, items);
				addPiece(piece); // TODO: removed true here, didnt compiled. Test.
			}
			setValueAtCurrent(f({ScalarType(t1)}));
		} catch (std::logic_error &e) {
			throw rethrow("DDEPiecewisePolynomialCurve::__construct__(grid,t0,t1,order,capd::map): ", e);
		}
	}

	/** number of pieces */
	size_type length() const { return m_pieces.size(); }

	/**
	 * returns subcurve from index1 to index2 (in indexaction of internal data. Mainly for internal use.
	 * Better use subcurve(TimePoint, TimePoint) version in your programs.
	 */
	Class subcurve(size_type index1, size_type index2) const {
		TimePointType from = (index1 == m_pieces.size() ? t0() : m_pieces[index1]->t0());
		Class result(m_grid, from, dimension());
		for (size_type i = index1; i < index2; ++i){
			result.addPiece(*m_pieces[i]);
		}
		result.setValueAtCurrent(index2 == m_pieces.size()? this->m_valueAtCurrent : (*m_pieces[index2])[0]);
		return result;
	}
	/**
	 * Returns subcurve from idex1 to current time. Mainly for internal use.
	 * Better use subcurve(TimePoint, TimePoint) version in your programs.
	 */
	Class subcurve(size_type index1) const { return subcurve(index1, m_pieces.size()); }

	/** returns subcurve of the segment between t0 and t1 */
	Class subcurve(TimePointType const& t0, TimePointType const& t1) const {
		size_type index1, index2;
		try { index1 = pointToIndex(t0); } catch (std::domain_error &e) { throw rethrow(badge() + "::subcurve(t0, t1): t0 outside range. ", e); }
		try { index2 = pointToIndex(t1); } catch (std::domain_error &e) { throw rethrow(badge() + "::subcurve(t0, t1): t1 outside range. ", e); }
		return subcurve(index1, index2);
	}

	/** returns subcurve of the segment between t0 and currentTime() */
	Class subcurve(TimePointType const& t0) const {
		size_type index1;
		try { index1 = pointToIndex(t0); } catch (std::domain_error &e) { throw rethrow(badge() + "::subcurve(t0): t0 outside range. ", e); }
		return subcurve(index1, m_pieces.size());
	}

	/** todo: docs */
	reverse_iterator rbegin() { return m_pieces.rbegin(); }
	/** todo: docs */
	const_reverse_iterator rbegin() const { return m_pieces.rbegin(); }
	/** todo: docs */
	reverse_iterator rend() { return m_pieces.rend(); }
	/** todo: docs */
	const_reverse_iterator rend() const { return m_pieces.rend(); }
	/** todo: docs */
	iterator begin() { return m_pieces.begin(); }
	/** todo: docs */
	const_iterator begin() const { return m_pieces.begin(); }
	/** todo: docs */
	iterator end() { return m_pieces.end(); }
	/** todo: docs */
	const_iterator end() const { return m_pieces.end(); }
	/** todo: docs */
	iterator at(TimePointType const & t) {
		// TODO: (RETHINK, FUTURE)  this is O(n) access, but deque does not allow to get hold on an iterator other than begin, end... :(
		iterator i = this->begin();
		while (i != this->end()){
			if ((**i).t0() == t) return i; else ++i;
		}
		std::ostringstream info;
		info << badge() << "::at(t) const: ";
		info << "time " << t << " not found. ";
		info << "Solution range is " << RealType(pastTime()) << " to " << RealType(t0());
		throw std::range_error(info.str());
		//		// below is bad code, due to the lack of access in deque
//		try { return m_pieces[pointToIndex(t)]; }
//		catch (std::domain_error &e) { throw rethrow("DDESolutionCurve::at(t): ", e); }
	}
	/** todo: docs */
	const_iterator at(TimePointType const & t) const {
		// TODO: (RETHINK, FUTURE) this is O(n) access, but deque does not allow to get hold on an iterator other than begin, end... :(
		const_iterator i = this->begin();
		while (i != this->end()){
			if ((**i).t0() == t) return i; else ++i;
		}
		std::ostringstream info;
		info << badge() << "::at(t) const: ";
		info << "time " << t << " not found. ";
		info << "Solution range is " << RealType(pastTime()) << " to " << RealType(t0());
		throw std::range_error(info.str());
//		// below is bad code, due to the lack of access in deque
//		try { return m_pieces[pointToIndex(t)]; }
//		catch (std::domain_error &e) { throw rethrow("DDESolutionCurve::at(t) const: ", e); }
	}

	/** todo: docs */
	GridType const& grid() const { return m_grid; }
	/** todo: docs */
	size_type dimension() const { return m_dimension; }
	/** todo: docs */
	size_type storageDimension() const {
		size_type dim = m_valueAtCurrent.dimension();
		for (auto it = begin(); it != end(); ++it){
			dim += (*it)->storageDimension();
		}
		return dim;
	}

	/** changes the time at which the curve ends on the right. WARNING: this is mostly for internal use, should not be changed directly */
	Class& setCurrentTime(TimePointType const& t0) { m_t_current = t0; return *this; }
	/** Returns the current time of the curve (right domain)  */
	TimePointType getCurrentTime() const { return m_t_current; }
	/** Return the time point of the first jet in the curve (left domain) */
	TimePointType getPastTime() const { return m_pieces.size() ? (**m_pieces.begin()).t0() : t0(); }

	/** same as getCurrentTime(), TODO: (NOT URGENT) deprecated? left for backward compatibility? */
	TimePointType t0() const { return getCurrentTime(); }
	/** same as getCurrentTime(), TODO: (NOT URGENT) deprecated? left for backward compatibility? */
	TimePointType getT0() const { return getCurrentTime(); }
	/** same as setCurrentTime(), TODO: (NOT URGENT) deprecated? left for backward compatibility? */
	Class& setT0(TimePointType const& t0) { return this->setCurrentTime(t0); }
	/** CAPD map compatible, same as getPastTime() */
	TimePointType leftDomain() const { return getPastTime(); } ///< CAPD map compatible
	/** CAPD map compatible, same as getCurrentTime() */
	TimePointType rightDomain() const { return t0(); } ///<
	/** TODO: (NOT URGENT) deprecated? left for backward compatibility? */
	TimePointType pastTime() const { return getPastTime(); }
	/** TODO: (NOT URGENT) deprecated? left for backward compatibility? */
	TimePointType currentTime() const { return getCurrentTime(); }

	/**
	 * returns raw representation of this solution as a vector.
	 * Use with caution, as it forgets all the structure of the solution.
	 * See @method get_x() for more information.
	 */
	operator VectorType() const { return this->get_x(); }

	/**
	 * returns raw representation of this solution as a vector.
	 * Use with caution, as it forgets all the structure of the solution.
	 *
	 * The representation is as follows:
	 *
	 * - first d coordinates: value at t=t0 = currentTime() (i.e. x(0) \in \R^d).
	 * - next, the jet at t = t0 - grid.h(), first the 0-th element of the jet, then order 1, and so on. Each element of the jest is d dimensional.
	 * - the following jets, up the the last piece.
	 * This might be a big vector, its size may vary, as in the representation, each jet can have in theory have different order.
     */
	VectorType get_x() const {
		VectorType x(storageDimension());
		size_type d = dimension();
		VectorType part = VectorType(m_valueAtCurrent);
		size_type I = 0;
		for (size_type i = 0; i < d; ++i, ++I)
			x[I] = part[i];
		for (auto ijet = rbegin(); ijet != rend(); ++ijet)
			for (auto icoeff = (*ijet)->begin(); icoeff != (*ijet)->end(); ++icoeff){
				part = VectorType(*icoeff);
				for (size_type i = 0; i < d; ++i, ++I)
					x[I] = part[i];
			}
		return x;
	}
	/** todo: docs */
	Class& set_x(VectorType const& x) {
		if (x.dimension() < storageDimension())
			throw std::logic_error("DDEPiecewisePolynomialCurve::set_x(): dimension of x not enough to fill all data in the curve.");
		if (x.dimension() > storageDimension())
			throw std::logic_error("DDEPiecewisePolynomialCurve::set_x(): dimension of x too big to fit data in the curve.");
		size_type d = dimension();
		VectorType part(d);
		size_type I = 0;
		for (size_type i = 0; i < d; ++i, ++I)
			part[i] = x[I];
		setValueAtCurrent(part);
		for (auto ijet = rbegin(); ijet != rend(); ++ijet)
			for (auto icoeff = (*ijet)->begin(); icoeff != (*ijet)->end(); ++icoeff){
				for (size_type i = 0; i < d; ++i, ++I)
					part[i] = x[I];
				(*icoeff) = part;
			}
		return *this;
	}

	/** todo: docs */
	Class& affineTransform(MatrixType const &M, VectorType const &v){
		throw std::logic_error(badge() + "::affineTransform(): Not Implemented Yet");
	}
	/** todo: docs */
	Class& translate(VectorType const &v){
		throw std::logic_error(badge() + "::translate(): Not Implemented Yet");
	}
	/** todo: docs */
	Class& add(VectorType const & v){
		throw std::logic_error(badge() + "::add(): Not Implemented Yet");
	}
	/** todo: docs */
	Class& add(Class const & set){
		throw std::logic_error(badge() + "::add(): Not Implemented Yet");
	}
	/** todo: docs */
	Class& mulThenAdd(ScalarType const& c, Class const & set){
		throw std::logic_error(badge() + "::mulThenAdd(): Not Implemented Yet");
	}
	/** todo: docs */
	Class& mul(ScalarType const& c){
		for (auto ijet = begin(); ijet != end(); ++ijet)
			for (auto icoeff = (*ijet)->begin(); icoeff != (*ijet)->end(); ++icoeff)
				(*icoeff) *= c;
		m_valueAtCurrent *= c;
		return *this;
	}
	/** todo: docs */
	Class& operator*=(ScalarType const& c){
		return this->mul(c);
	}

	/** returns index in the array of the Jets for a given TimePoint */
	size_type pointToIndex(TimePointType const& t) const {
		int ibase = int(pastTime());
		int i0 = int(currentTime()) - ibase;
		int i = int(t) - ibase;
		if (i < 0) 		throw std::domain_error(badge() + "::pointToIndex(): piece index before the past jet.");
		if (i >= i0) 	throw std::domain_error(badge() + "::pointToIndex(): piece index bigger than number of jets.");
		return i;
	}
	CurvePieceType& getPiece(TimePointType const& t) const {
		size_type i;
		try { i = pointToIndex(t); } catch (std::domain_error &e) { throw rethrow(badge() + "::getPiece(TimePoint):", e); }
		return *(m_pieces[i]);
	}
	/** it assumes that index is already processed (see pointToIndex()) */
	CurvePieceType& getPiece(size_type const& index) const {
		// TODO: (IMPORTANT, RETHINK): this will slow computations down - make it to not check the constraints or add compiler flag to disabkle extra checking for speed purposes.
		// TODO: (IMPORTANT, CHECK): check other time-critical code if it uses getPiece(TimePoint) or this one and make sure it uses this one, without checking.
		int ibase = int(pastTime());
		int i0 = int(currentTime()) - ibase;
		if (index < 0 || index >= i0)
			throw std::domain_error(badge() + "::getPiece(size_type): piece index outside range.");
		return *(m_pieces[index]);
	}

	// TODO: only implement eval for RealType (?) and make it rigorous in cases.
	VectorType eval(TimePointType t) const {
		try{
			return t == t0() ? VectorType(m_valueAtCurrent) : getPiece(t).eval(t);
		} catch (std::domain_error &e){
			throw rethrow(badge() + "::eval(TimePoint):", e);
		}
	}
	/**
	 * WARNING: THIS FUNCTION IS NOT RIGOROUS FOR INTERVALS! DO NOT USE IN PROOFS
	 * WARNING: COMPUTE RIGOROUSLY ONLY USING eval(TimePointType t) VERSION
	 */
	VectorType eval(RealType t) const {
		// TODO: (FUTURE???) raise warning (compilation time?) when this function is used?
		// TODO: (FUTURE/IMPORTANT?) it should be possible to make it rigorous, you need
		RealType epsi;
		TimePointType ti = m_grid.point(0);
		m_grid.split(t, ti, epsi);
		try{
			return (ti == t0() && epsi == 0.0) ? VectorType(m_valueAtCurrent) : getPiece(ti).evalAtDelta(epsi);
		} catch (std::domain_error &e){
			throw rethrow(badge() + "::eval(Real):", e);
		}
	}
	void eval(TimePointType t, DataType& out) const {
		try{
			out = (t == t0() ? m_valueAtCurrent : getPiece(t)[0]);
		} catch (std::domain_error &e){
			throw rethrow(badge() + "::eval(TimePoint, SetType&):", e);
		}
	}
	void eval(RealType t, DataType& out) const {
		throw std::logic_error(badge() + "::eval(RealType, SetType&): Not implemented yet.");
	}

	Class dt(DataType const& valueAtCurrent, size_type n = 1) const {
		if (n == 0) { return *this; } // no need to compute anything
		Class result(m_grid, pastTime(), dimension());
		for (auto ijet = begin(); ijet != end(); ++ijet)
			result.addPiece((*ijet)->dt(n));
		result.setValueAtCurrent(valueAtCurrent);
		return result;
	}

	template<typename FunctionalSpec>
	Class dt(FunctionalSpec const& f, size_type n = 1) const {
		if (n == 0) { return *this; } // no need to compute anything
		Class result(m_grid, pastTime(), dimension());
		for (auto ijet = begin(); ijet != end(); ++ijet)
			result.addPiece((*ijet)->dt(n));
		if (n == 1){
			// faster and simplified than the case n > 1
			DataType v = f(*this);
			result.setValueAtCurrent(v);
		}else{
			typedef typename FunctionalSpec::ValueStorageType ValueStorageType;
			ValueStorageType coeffs;
			f.computeDDECoefficients(*this, coeffs);
			DataType v = coeffs[n-1];
			while (n > 0) v *= n--; // v = x^(n)/n!, so we need to multiply by n!
			result.setValueAtCurrent(v);
		}
		return result;
	}

	Class increasedOrder(size_type r = 1){
		if (r == 0) { return *this; } // no need to compute anything
		Class result(m_grid, pastTime(), dimension());
		for (auto ijet = begin(); ijet != end(); ++ijet)
			result.addPiece((*ijet)->increasedOrder(r));
		result.setValueAtCurrent(m_valueAtCurrent);
		return result;
	}
	Class decreasedOrder(size_type r = 1){
		if (r == 0) { return *this; } // no need to compute anything
		Class result(m_grid, pastTime(), dimension());
		for (auto ijet = begin(); ijet != end(); ++ijet)
			result.addPiece((*ijet)->decreasedOrder(r));
		result.setValueAtCurrent(m_valueAtCurrent);
		return result;
	}

	/*
	 * add piece of curve by making a copy,
	 * the copy will be owned by the DDESolutionCurve
	 * returns this for a cascade
	 * THIS SHOULD BE PREFFERED WAY TO ADD PIECES BY END USER
	 */
	Class& addPiece(CurvePieceType const & newPiece) { addPiece(new CurvePieceType(newPiece)); return *this; }
	/** adds directly by pointer, faster than by copying from reference. The curve will be resp. for deleting. */
	Class& addPiece(CurvePieceType* newPiece) { newPiece->setT0(t0()); m_pieces.push_back(newPiece); setT0(t0() + getStep()); return *this; }
	/** same as addPiece but adds in the past (before existing pieces) */
	Class& addPastPiece(const CurvePieceType& newPiece) { addPastPiece(new CurvePieceType(newPiece)); return *this; }
	/** adds directly by pointer, faster than by copying from reference. The curve will be resp. for deleting. */
	Class& addPastPiece(CurvePieceType* newPiece) { newPiece->setT0(pastTime() - getStep()); m_pieces.push_front(newPiece); return *this; }
	/** set the value at a current time */
	Class& setValueAtCurrent(const DataType& value) { m_valueAtCurrent = value; return *this; }
	/** get the value at current */
	DataType& getValueAtCurrent() { return m_valueAtCurrent; };
	/** get the value at current (const object) */
	const DataType& getValueAtCurrent() const { return m_valueAtCurrent; };

	/** jet at grid point t0 and valid until t1 (copy, see j() for reference) */
	JetType jet(RealType t0, RealType t1) const { throw std::logic_error(badge() + "::jet: only supports jets over grid intervals."); }
	/** jet at grid point t0 and valid until t1 (copy, see j() for reference) */
	JetType jet(TimePointType t0, TimePointType t1) const {
		if (t0 + 1 != t1) throw std::logic_error(badge() + "::jet: does not support intervals longer than grid step size");
		try { return getPiece(t0); } catch (std::domain_error &e){ throw rethrow(badge() + "::jet(TimePoints):", e); }
	}
	/** jet at grid point t0, and valid for t0 to t0 + grid_step (copy, see j() for reference) */
	JetType jet(TimePointType t0) const {
		try { return getPiece(t0); } catch (std::domain_error &e){ throw rethrow(badge() + "::jet(TimePoint):", e); }
	}
	/** use with caution! Mainly for internal use for speed purposes */
	JetType* jetPtr(TimePointType t0) const {
		try { return &getPiece(t0); } catch (std::domain_error &e){ throw rethrow(badge() + "::jetPtr(TimePoint):", e); }
	}
	/** returns value of the curve at the time point i (size_type) */
	DataType const& value(size_type i) const { try { return m_grid(i) == m_t_current ? m_valueAtCurrent : getPiece(i)[0]; } catch (std::domain_error &e){ throw rethrow(badge() + "::value(index): ", e); } }
	/** returns value of the curve at the time point t0 (TimePointType) */
	DataType const& value(TimePointType t0) const { try { return t0 == m_t_current ? m_valueAtCurrent : getPiece(t0)[0]; } catch (std::domain_error &e){ throw rethrow(badge() + "::value(TimePoint): ", e); } }
	/** returns value of the curve at the time point i (size_type) */
	DataType& value(size_type i) { try { return m_grid(i) == m_t_current ? m_valueAtCurrent : getPiece(i)[0]; } catch (std::domain_error &e){ throw rethrow(badge() + "::value(index): ", e); } }
	/** returns value of the curve at the time point t0 (TimePointType) */
	DataType& value(TimePointType t0) { try { return t0 == m_t_current ? m_valueAtCurrent : getPiece(t0)[0]; } catch (std::domain_error &e){ throw rethrow(badge() + "::value(TimePoint): ", e); } }

	/** k-th coeff at grid point i */
	DataType const& j(size_type i, size_type k) const { try { return getPiece(i)[k]; } catch (std::domain_error &e){ throw rethrow(badge() + "::j(index, order): ", e); } }
	/** k-th coeff at grid point t0 */
	DataType const& j(TimePointType t0, size_type k) const { try { return getPiece(t0)[k]; } catch (std::domain_error &e){ throw rethrow(badge() + "::j(TimePoint, order): ", e); } }
	/** k-th coeff at grid point i */
	DataType& j(size_type i, size_type k) { try { return getPiece(i)[k]; } catch (std::domain_error &e){ throw rethrow(badge() + "::j(index, order): ", e); } }
	/** k-th coeff at grid point t0 */
	DataType& j(TimePointType t0, size_type k) { try { return getPiece(t0)[k]; } catch (std::domain_error &e){ throw rethrow(badge() + "::j(TimePoint, order): ", e); } }
	/** jet at grid point t0 (reference) */
	JetType const& j(TimePointType t0) const { try { return getPiece(t0); } catch (std::domain_error &e){ throw rethrow(badge() + "::j(TimePoint): ", e); } }
	/** jet at grid point t0 (reference) */
	JetType& j(TimePointType t0) { try { return getPiece(t0); } catch (std::domain_error &e){ throw rethrow(badge() + "::j(TimePoint): ", e); } }
	/** returns order of jet available at point t. This is for optimization purposes (do not copy whole Jet just to know its order) */
	size_type jetOrderAt(TimePointType t0){
		try { return getPiece(t0).order(); } catch (std::domain_error &e){ throw rethrow(badge() + "::jetOrderAt(TimePoint): ", e); }
	}
	/** returns maximal common jet of all jets in this curve  */
	size_type jetCommonMaximalOrder(){
		try {
			// TODO: use std::algorithm to do that DRY?
			auto it = m_pieces.begin();
			size_type min_order = (*it)->order(); ++it;
			for (; it != m_pieces.end(); ++it){
				auto order = (*it)->order();
				min_order = (min_order > order ? order : min_order);
			}
			return min_order;
		} catch (std::domain_error &e){
			throw rethrow(badge() + "::jetCommonMaximalOrder(): ", e);
		}
	}

	/** returns grid object for this curve */
	const GridType& getGrid() const { return m_grid; }
	/** returns the grid step for this curve */
	const TimePointType getStep() const { return m_grid.point(1); }

	/** show human readable representation */
	std::string show() const;
	/** friend operator must be inline (in fact it does not, but then I have to move this to .cpp to avoid linker error... I leave it as it is for now. TODO:), or the linker gets confused, so we use writeTo/readFrom */
	friend std::ostream& operator<<(std::ostream& out, Class const& curve){ curve.writeTo(out); return out; }
	/** friend operator must be inline (in fact it does not, but then I have to move this to .cpp to avoid linker error... I leave it as it is for now. TODO:), or the linker gets confused, so we use writeTo/readFrom */
	friend std::istream& operator>>(std::istream& in, Class const& curve){ curve.readFrom(in); return in; }

	/**
	 * Naming convention. To be a pair of extend <-> epsilonShift
	 * See move() for more docs.
	 */
	template<typename DynSysSpec>
	Class& extend(DynSysSpec& solver) {
		throw std::logic_error(badge() + "::extend(): Not Implemented Yet");
	}

	/**
	 * it computes epsilon shift of in_curve and store it in out_result.
	 * out_result must be setup so that all desired TimePoints
	 * are presented and Jets must be initialized to a desired order.
	 * Most probable scenario: out_result has Jets of order n
	 * at times -tau, -tau+1/h, ..., -h, 0 (current time), where h =
	 * 1/p, as in the paper. Then we shift the final part of the in_curve
	 * starting at t0 - tau - h up to -h, to t0 - tau - h + epsi, ..., -h + epsi.
	 * The coefficients will be stored in out_result at -tau, .., 0, respectively!
	 * The t0 is taken from in_curve.t0() in this procedure. For a general t0
	 * see other function.
	 *
	 * TODO: (FUTURE) this should be not templated by DynSysSpec but concrete "interfaced" type, known in advance (?)
	 */
	template<typename DynSysSpec>
	void epsilonShift(
				DynSysSpec const& solver,
				RealType const& epsilon,
				Class& out_result) const {
		epsilonShift(solver, t0() - 1, epsilon, out_result);
	}

	/**
	 * it computes epsilon shift of in_curve and store it in out_result.
	 * out_result must be setup so that all desired TimePoints
	 * are presented and Jets must be initialized to a desired order.
	 * Most probable scenario: out_result has Jets of order n
	 * at times -tau, -tau+1/h, ..., -h, 0 (current time), where h =
	 * 1/p, as in the paper. Then we shift the final part of the in_curve
	 * starting at at_t0 - tau - h up to -h, to t0 - tau - h + epsi, ..., -h + epsi.
	 * The coefficients will be stored in out_result at -tau, .., 0, respectively!
	 *
	 * NOTE: epsilon should be >= 0!
	 *
	 * TODO: (FUTURE) this should be not templated by DynSysSpec but concrete "interfaced" type, known in advance
	 */
	template<typename DynSysSpec>
	void epsilonShift(
				DynSysSpec const& solver,
				TimePointType const& at_t0, RealType const& epsilon,
				Class& out_result) const {
		// this is only forward epsilon step. I have tried implementing backward step,
		// but it was a disaster, due to the representation mismatch at grid points.
		// backward epsi step might be better when a smooth representation is available!
		// TODO: (FUTURE) use solver and history to recompute Xi over shorter intervals t_i + [0, \epsi], when computing coefficients?
		// TODO: (NOT URGENT) add checks for sanity of the input before proceeding.
		TimePointType how_far = out_result.t0() - out_result.pastTime();
		auto ijet_out = out_result.begin();
		auto ijet_this = this->at(at_t0 - how_far);
		for (; ijet_out != out_result.end(); ++ijet_out, ++ijet_this){
			(**ijet_this).evalAtDelta(epsilon, (**ijet_out)[0]);
			size_type order = (**ijet_out).order();
			for (size_type k = 1; k <= order; ++k){
				(**ijet_this).evalCoeffAtDelta(k, epsilon, (**ijet_out)[k]);
			}
		}
		// we assume this extra jet is already computed,
		// otherwise, we would need to call solver to extend solution - we do not want this.
		(**ijet_this).evalAtDelta(epsilon, out_result.m_valueAtCurrent);
	}

	/** computes dot taking appropriate number of first components from Curve */
	ScalarType dot(VectorType const& v) const{
		throw std::logic_error(badge() + "::dot(): Not Implemented Yet");
	}
	/** computes dot taking appropriate jets into account. TODO: only template with a given type of section (based on JetType)? */
	template<typename JetSection>
	ScalarType dot(JetSection const& v) const {
		auto ijet = this->rbegin();
		ScalarType result = VectorType(this->m_valueAtCurrent) * v.get_s();
		for (auto iv = v.rbegin(); iv != v.rend(); ++iv, ++ijet){
			auto icoeff = (*ijet)->begin();
			for (auto ivc = iv->begin(); ivc != iv->end(); ++ivc, ++icoeff){
				result += VectorType(*icoeff) * VectorType(*ivc);
			}
		}
		return result;
	}

	/** clears the Curve, removing all pieces and setting it to be a point of value 0. at current t0 */
	Class& clear(){
		m_valueAtCurrent *= 0.;
		deallocatePieces();
		return *this;
	}

	/** standard C++ thing */
	virtual ~DDEPiecewisePolynomialCurve() { this->deallocatePieces(); }
protected:
	const GridType& m_grid;
	TimePointType m_t_current;
	VectorType m_lastEnclosure;
	DataType m_valueAtCurrent;
	size_type m_dimension;
	PiecesStorageType m_pieces;

	/** safely deallocates pointers in m_pieces and clears the container */
	void deallocatePieces(){
		auto pIt = begin();
		for (; pIt != end(); ++pIt){
			if (*pIt) delete *pIt;
			// delete *pIt is ok here, as *pIt is a CurvePieceTypePtr
			// (**pIt is the CurvePiece object, sic!).
			*pIt = NULL;
		}
		m_pieces.clear();
	}

	/** low level, does not care about changing time points, does not clear m_pieces, etc. just copies from to */
	template<typename IteratorSpec>
	void copyPieces(IteratorSpec from, IteratorSpec to){
		// TODO: (FUTURE) that's not the most optimal implementation...
		for (; from != to; ++from)
			m_pieces.push_back(new CurvePieceType(*from));
	}
	/** low level, does not care about changing time points, does not clear m_pieces, etc. just copies from to */
	void copyPieces(iterator from, iterator to){
		// TODO: (NOT URGENT) that's not the most optimal implementation...
		for (; from != to; ++from)
			m_pieces.push_back(new CurvePieceType(**from));
	}
	/** low level, does not care about changing time points, does not clear m_pieces, etc. just copies from to */
	void copyPieces(const_iterator from, const_iterator to){
		// TODO: (NOT URGENT) that's not the most optimal implementation...
		for (; from != to; ++from)
			m_pieces.push_back(new CurvePieceType(**from));
	}

	/** internal function to make short impl. of << operator and the long impl. in .hpp file (see comment in << operator) */
	void writeTo(std::ostream& out) const;
	/** internal function to make short impl. of >> operator and the long impl. in .hpp file (see comment in >> operator) */
	void readFrom(std::istream& out) const;
};


} // namespace ddes
} // namespace capd


#endif /* _CAPD_DDES_DDEPIECEWISEPOLYNOMIALCURVE_H_ */
