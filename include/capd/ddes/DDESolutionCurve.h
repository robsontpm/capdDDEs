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

#ifndef _CAPD_DDES_DDESOLUTIONCURVE_H_
#define _CAPD_DDES_DDESOLUTIONCURVE_H_

#include <capd/basicalg/TypeTraits.h>
#include <stdexcept>
#include <string>
#include <deque>
#include <capd/ddes/DDECommon.h>
#include <capd/ddes/DDEForwardTaylorCurvePiece.h>

namespace capd{
namespace ddes{

/**
 * A class to represent the solution $x(t)$ of the DDE. It provides methods
 * To ask for the DDECurvePiece over specific intervals (time intervals)
 * The basic implementation allow to ask about intervals $I_i = [i * h, i* h + h)$,
 * where $h = {\tau}\over{p}$ is the grid size as in FoCM article / PhD thesis.
 * The DDECurvePiece is just a Jet of order n of this solution x at a given time
 * and valid on the specified interval + estimates on the n+1-st derivative.
 * Please consult DDE(ForwardTaylor)CurvePiece for more information how
 * this is organized.
 *
 * DDESolutionCurve also contains the value of the solution at a maximal time
 * $t_{current}$ - to be used later in the process of solving DDEs.
 * Usually, we will start with some candidate DDESolutionCurve defined over
 * [-\tau, 0], with $t_{current} = 0$. Then we propagate solution
 * by Integrator (see DDESolver) to some $t_{future}$, so that
 * the solution is estimated over [-\tau, t_{future}]. Please
 * note that the history is preserved for theoretical and technical purposes.
 *
 * From the theoretical point of view, the history could be used to design
 * solvers for varying delay (as long as it is bounded).
 * From technical point of view, there might be some ways to improve current
 * methods using the history (out of scope as of now).
 *
 * TODO: (NOT URGENT, RETHINK) see the comments about allowing change of m_grid in DDEPiecewisePolynomialCurve
 * TODO: (DRY, RETHINK) now a lot of code is duplicated from DDEPIecewisePolynomiaCurve, this class should inherit from there maybe and implement extra! (The problem here is that, we have extra Lohner set structure here, and how to make it compatible with DDEPIecewisePolynomial which works on simple vectors?)
 * TODO: (FAR FUTURE) see notes below this line
 * TODOs / DevNOTEs:
 *
 * Newer note:
 * As of now, class is designed to work with my own implementation of Lohner sets
 * (see 'storage/' subdirectory). They are just Doubletons with some restrictions
 * of classical CAPD sets removed to allow m_C and m_r0 parts to be more freely
 * defined. Ideally I would like to use standard sets from CAPD, but this
 * might be not possible / not optimal (because of the part r_0 in definitions,
 * which is ideally shared between sets). Alternatively, I might consider
 * using a big CAPD set and some proxies that can extract subcolumns/subrows
 * of various elements to be passed down to ForwardTaylorCurvePieces or
 * Taylor coefficients of given order.
 *
 * Older note:
 * The class is templated with a classical CAPD Lohner Set Type, but we will
 * forcefully alter it behavior (I hope). I will probably use its QR policy
 * etc., but I will alter the matrix m_C and vector m_r - they will
 * be stored from within the set more or less, and might not be related to
 * the space dimension. What is required that the product $m_C * m_r$ has
 * dimension $d$, i.e. m_C is $d$ rows, $N0$ columns, and $m_r0$ is of dimension $N0$.
 */
template<typename SetSpec>
class DDESolutionCurve :
		public DoubletonInterface<typename SetSpec::MatrixType, typename SetSpec::Policy>{
public:

	typedef SetSpec SetType;
	typedef SetType DataType;
	typedef DoubletonInterface<typename SetSpec::MatrixType, typename SetSpec::Policy> BaseClass;
	typedef DDESolutionCurve<SetSpec> Class;
	typedef typename SetType::MatrixType MatrixType;
	typedef typename MatrixType::RowVectorType VectorType;
	typedef typename MatrixType::ScalarType ScalarType;
	typedef ScalarType RealType; // TODO: (FAR FUTURE, NOT URGENT) do as in capd by traits
	// TODO: (FAR FUTURE) allow to change the grid type
	typedef DiscreteTimeGrid<RealType> GridType;
	typedef typename GridType::TimePointType TimePointType;
	typedef typename BaseClass::size_type size_type;
	// TODO: (FAR FUTURE) allow to change the curve piece type
	typedef DDEForwardTaylorCurvePiece<TimePointType, SetType> CurvePieceType;
	typedef typename SetType::Policy QRPolicy;  // I prefer QRPolicy name....
	typedef typename SetType::Policy Policy;    // to match CAPD naming convention
	typedef CurvePieceType JetType;
	typedef std::deque<CurvePieceType*> PiecesStorageType;
	typedef std::deque<bool> OwnerStorageType;
	typedef typename PiecesStorageType::iterator iterator;
	typedef typename PiecesStorageType::const_iterator const_iterator;
	typedef typename PiecesStorageType::reverse_iterator reverse_iterator;
	typedef typename PiecesStorageType::const_reverse_iterator const_reverse_iterator;

	/** copy constructor, standard thing. */
	DDESolutionCurve(const DDESolutionCurve& orig):
			m_grid(orig.m_grid), m_t_current(orig.m_t_current),
			m_r0(new VectorType(*orig.m_r0)),
			m_valueAtCurrent(orig.m_valueAtCurrent),
			m_dimension(orig.m_dimension)
	{
		copyPieces(orig.begin(), orig.end());
		updateCommonR0();
	}

	/** assign operator. NOTE: it wont work if the @param orig has different grid! It will throw std::logic_error(). */
	DDESolutionCurve& operator=(const DDESolutionCurve& orig){
		if (&m_grid != &(orig.m_grid))
			throw std::logic_error("DDESolutionCurve::operator=(orig): cannot assign solution over a different grid.");
		try {
			deallocate();
			allocate(orig.storageN0());
			*m_r0 = *(orig.m_r0);
			m_t_current = orig.m_t_current;
			m_dimension = orig.m_dimension;
			m_valueAtCurrent = orig.m_valueAtCurrent;
			copyPieces(orig.begin(), orig.end());
			updateCommonR0();
		} catch (std::logic_error& e){
			throw rethrow("DDESolutionCurve::operator=(orig): ", e);
		}
		return *this;
	}

	/** makes a Curve that is d dimensional, over given grid, located initially at TimePoint 0, and its r0 vector is N0 dimensional */
	DDESolutionCurve(const GridType& grid, size_type d = 0, size_type N0 = 0):
			m_grid(grid), m_t_current(grid.point(0)),
			m_r0(new VectorType(N0)),
			m_valueAtCurrent(VectorType(d), m_r0),
			m_dimension(d)
	{ /* no need to update, no pieces, current updated in constructor */ }

	/** makes a Curve that is d dimensional, over given grid starting at a given TimePoint, and its r0 vector is N0 dimensional */
	DDESolutionCurve(const TimePointType& t0, size_type d = 0, size_type N0 = 0):
			m_grid(t0.grid()), m_t_current(t0),
			m_r0(new VectorType(N0)),
			m_valueAtCurrent(VectorType(d), m_r0),
			m_dimension(d)
	{ /* no need to update, no pieces, current updated in constructor */ }

	/** creates a Curve that is a single point of given value at time point t0,m the r0 part of value is used as r0 of the Curve */
	DDESolutionCurve(const TimePointType& t0, const SetType& value):
			m_grid(t0.grid()), m_t_current(t0),
			m_r0(new VectorType(value.get_r0())),	// I do not assume in DoubletonInterface, that I can transfer ownership (Yet, TODO: (NOT URGENT) rething take_*() semantics and implement it simply in BasicDoubleton as return new Vector(m_r0), it should work)
			m_valueAtCurrent(value),				// copy set structure. we will need to update r0 (se above)
			m_dimension(VectorType(value).dimension())
	{ try { updateCommonR0(); } catch (std::logic_error& e) { throw rethrow("DDESolutionCurve::__construct__(grid,t0,value): ", e); } }

	/**
	 * creates a Curve that is a single point of given value at time point t0 and with given dimension of r0 part.
	 * It splits x = value vector into the set mid(x) + C * r0, C = Id, r0 = x - mid(x), as in CAPD.
	 * It will use this r0 for the Lohner structure for the subsequent jets, i.e. when doing ODE integration.
	 */
	DDESolutionCurve(const TimePointType& t0, const VectorType& value):
			m_grid(t0.grid()), m_t_current(t0),
			m_r0(0),
			m_valueAtCurrent(value),
			m_dimension(value.dimension())
	{
		m_r0 = m_valueAtCurrent.take_r0();
		// no need to update, no pieces, current updated in constructor
	}

	/**
	 * creates a constant function with a given value on interval [t1, t0]
	 * and with given order (useful for future use as a a set in integration).
	 * The r0 part will be COPIED from the set and will be common to all elements.
	 */
	DDESolutionCurve(
				const TimePointType& t0,
				const TimePointType& t1,
				size_type order,
				const SetType& value):
			m_grid(t0.grid()), m_t_current(t0),
			m_r0(new VectorType(value.get_r0())),
			m_valueAtCurrent(value, m_r0),
			m_dimension(VectorType(value).dimension())
	{
		for (TimePointType t = t0; t < t1; ++t)
			addPiece(new CurvePieceType(t, order, m_valueAtCurrent, m_r0), true); // is faster than by reference, true => pass the ownership
		// do not need to updateCommonR0(), we used constructors to assure it.
	}

	/**
	 * creates a constant function with a given value on interval [t1, t0]
	 * and with given order (useful for future use as a a set in integration).
	 * The r0 part will be transferred (take from, the pointer - of possible / the set implementation allows that)
	 * from the set and will be common to all elements.
	 */
	DDESolutionCurve(
				const TimePointType& t0,
				const TimePointType& t1,
				size_type order,
				SetType& value):
			m_grid(t0.grid()), m_t_current(t0),
			m_r0(value.take_r0()),
			m_valueAtCurrent(value),
			m_dimension(VectorType(value).dimension())
	{
		try {
			updateCommonR0(); // in m_valueAtCurrent...
		}catch (std::logic_error& e) {
			throw rethrow("DDESolutionCurve::__construct__(grid,t0,t1,order,set): ", e);
		}
		for (TimePointType t = t0; t < t1; ++t)
			addPiece(new CurvePieceType(t, order, m_valueAtCurrent, m_r0), true); // is faster than by reference, true => pass the ownership
		// do not need to updateCommonR0(), we used constructors to assure it.
	}

	/**
	 * creates a constant function with a given value on interval [t1, t0]
	 * and with given order (useful for future use as a a set in integration).
	 * The r0 part will be N0 dimensional (default = 0, i.e. no C*r0 set part).
	 */
	DDESolutionCurve(
				const TimePointType& t0,
				const TimePointType& t1,
				size_type order,
				const VectorType& value,
				size_type N0 = 0):
			m_grid(t0.grid()), m_t_current(t0),
			m_r0(new VectorType(N0)),
			m_valueAtCurrent(value, m_r0),
			m_dimension(value.dimension())
	{
		for (TimePointType t = t0; t < t1; ++t)
			addPiece(new CurvePieceType(t, order, value, m_r0), true); // is faster than by reference, true => pass the ownership

		try {
			updateCommonR0();
		} catch (std::logic_error& e) {
			throw rethrow("DDESolutionCurve::__construct__(grid,t0,t1,order,vector,N0): ", e);
		}
	}
	/**
	 * creates a representation of the function f over the interval [t0,t1]
	 * The r0 part will be N0 dimensional (default = 0, i.e. no C*r0 set part).
	 *
	 * TODO: (RETHINK) only support for the same MatrixType as the internal MatrixType? Better for static compilation. Might include FMapType in the list of datattypes of this curve.... On the other hand, this way we might mix IMap with our custom Solution<MyMatrix>...
	 * TODO: (DRY, RETHINK) now we have some overlap with DDEPiecewisePolynomial (Nonrigorous version), cleanup the code and make the good inheritance structure?
	 */
	template<typename AnyMatrixSpec>
	DDESolutionCurve(
				const TimePointType& t0,
				const TimePointType& t1,
				size_type order,
				const capd::map::Map<AnyMatrixSpec>& f,
				size_type N0 = 0):
			m_grid(t0.grid()), m_t_current(t0),
			m_r0(new VectorType(N0)),
			m_valueAtCurrent(VectorType(f.imageDimension()), m_r0),
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
		if (f.degree() < order+2){
			std::ostringstream info;
			info << "DDESolutionCurve::__construct__(): f should have degree at least " << order+2 << ", ";
			info << "now f is of order " << f.degree() << ". ";
			info << "Consider using f.setOrder(order+2) on the parameter f, before calling this constructor. ";
			throw std::logic_error(info.str());
		}
		try {
			for (TimePointType t = t0; t < t1; ++t){
				FJet tt(f.imageDimension(), 1, order + 2);
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
				ScalarType HH(0., ScalarType(m_grid.point(1)).rightBound()); // TODO: this is kind of bad code...
				tti = tt.begin(0);
				*(tti) = ScalarType(t) + HH;
				*(++tti) = 1.;
				FJet FtXi = f(tt);
				VectorType xi(m_dimension);
				for (size_type i = 0; i < f.imageDimension(); ++i){
					auto ftxi = FtXi.begin(i) + order + 1;
					xi[i] = *ftxi;
				}
				CurvePieceType* piece = new CurvePieceType(t, items, N0);
				addPiece(piece, true); // is faster than by reference, true => pass the ownership
			}
			setValueAtCurrent(f({ScalarType(t1)}));
			MatrixType C(storageDimension(), storageN0());
			VectorType r0(storageN0());
			this->set_Cr0(C, r0);
			updateCommonR0();
		} catch (std::logic_error &e) {
			throw rethrow("DDESolutionCurve::__construct__(grid,t0,t1,order,capd::map,N0): ", e);
		}
	}

	/** number of pieces */
	size_type length() const { return m_pieces.size(); }

	/**
	 * returns subcurve from index1 to index2 (in indexaction of internal data. Mainly for internal use.
	 * Better use subcurve(TimePoint, TimePoint) version in your programs.
	 */
	DDESolutionCurve subcurve(size_type index1, size_type index2) const {
		TimePointType from = (index1 == m_pieces.size() ? t0() : m_pieces[index1]->t0());
		DDESolutionCurve result(from, dimension(), storageN0());
		result.set_r0(get_r0());
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
	DDESolutionCurve subcurve(size_type index1) const { return subcurve(index1, m_pieces.size()); }

	/** returns subcurve of the segment between t0 and t1 */
	DDESolutionCurve subcurve(TimePointType const& t0, TimePointType const& t1) const {
		size_type index1, index2;
		try { index1 = pointToIndex(t0); } catch (std::domain_error &e) { throw rethrow("DDESolutionCurve::subcurve(t0, t1): t0 outside range. ", e); }
		try { index2 = pointToIndex(t1); } catch (std::domain_error &e) { throw rethrow("DDESolutionCurve::subcurve(t0, t1): t1 outside range. ", e); }
		return subcurve(index1, index2);
	}

	/** returns subcurve of the segment between t0 and currentTime() */
	DDESolutionCurve subcurve(TimePointType const& t0) const {
		size_type index1;
		try { index1 = pointToIndex(t0); } catch (std::domain_error &e) { throw rethrow("DDESolutionCurve::subcurve(t0): t0 outside range. ", e); }
		return subcurve(index1, m_pieces.size());
	}

	reverse_iterator rbegin() { return m_pieces.rbegin(); }
	const_reverse_iterator rbegin() const { return m_pieces.rbegin(); }
	reverse_iterator rend() { return m_pieces.rend(); }
	const_reverse_iterator rend() const { return m_pieces.rend(); }
	iterator begin() { return m_pieces.begin(); }
	const_iterator begin() const { return m_pieces.begin(); }
	iterator end() { return m_pieces.end(); }
	const_iterator end() const { return m_pieces.end(); }
	iterator at(TimePointType const & t) {
		// TODO: (RETHINK, FUTURE)  this is O(n) access, but deque does not allow to get hold on an iterator other than begin, end... :(
		iterator i = this->begin();
		while (i != this->end()){
			if ((**i).t0() == t) return i; else ++i;
		}
		std::ostringstream info;
		info << "DDESolutionCurve::at(t) const: ";
		info << "time " << t << " not found. ";
		info << "Solution range is " << RealType(pastTime()) << " to " << RealType(t0());
		throw std::range_error(info.str());
		//		// below is bad code, due to the lack of access in deque
//		try { return m_pieces[pointToIndex(t)]; }
//		catch (std::domain_error &e) { throw rethrow("DDESolutionCurve::at(t): ", e); }
	}
	const_iterator at(TimePointType const & t) const {
		// TODO: (RETHINK, FUTURE) this is O(n) access, but deque does not allow to get hold on an iterator other than begin, end... :(
		const_iterator i = this->begin();
		while (i != this->end()){
			if ((**i).t0() == t) return i; else ++i;
		}
		std::ostringstream info;
		info << "DDESolutionCurve::at(t) const: ";
		info << "time " << t << " not found. ";
		info << "Solution range is " << RealType(pastTime()) << " to " << RealType(t0());
		throw std::range_error(info.str());
//		// below is bad code, due to the lack of access in deque
//		try { return m_pieces[pointToIndex(t)]; }
//		catch (std::domain_error &e) { throw rethrow("DDESolutionCurve::at(t) const: ", e); }
	}

	GridType const& grid() const { return m_grid; }
	size_type dimension() const { return m_dimension; }
	TimePointType t0() const { return m_t_current; }
	TimePointType currentTime() const { return t0(); }
	TimePointType rightDomain() const { return t0(); } ///< CAPD map compatible
	TimePointType pastTime() const { return m_pieces.size() ? (**m_pieces.begin()).t0() : t0(); }
	TimePointType leftDomain() const { return pastTime(); } ///< CAPD map compatible
	/** DoubletonStorageInterface: returns a dimension of the Jet as a Vector (sequence) to store all the coefficients (without Xi part) */
	size_type storageDimension() const {
		// size_type dim = m_valueAtCurrent.storageDimension();
		size_type dim = m_valueAtCurrent.dimension();
		for (auto it = begin(); it != end(); ++it){
			dim += (*it)->storageDimension();
		}
		return dim;
	}
	/** DoubletonStorageInterface: returns a dimension of the r0 vector */
	size_type storageN0() const { return m_r0->dimension(); }
	/** see DoubletonInterface */
	BaseClass& affineTransform(MatrixType const &M, VectorType const &v);
	/** see DoubletonInterface */
	BaseClass& translate(VectorType const &v);
	/** see DoubletonInterface */
	VectorType get_x() const;
	/** see DoubletonInterface */
	MatrixType get_C() const;
	/** see DoubletonInterface */
	VectorType get_r0() const;
	/** see DoubletonInterface */
	MatrixType get_B() const;
	/** see DoubletonInterface */
	VectorType get_r() const;
	/** see DoubletonInterface */
	MatrixType get_Binv() const;
	/** see DoubletonInterface */
	BaseClass& set_x(VectorType const &x);
	/** see DoubletonInterface */
	BaseClass& set_C(MatrixType const &C);
	/** see DoubletonInterface */
	BaseClass& set_r0(VectorType const &r0);
	/**
	 * Set C * r0 part of this set (as a doubleton). Because this set is defined by many sub-sets
	 * then we need to assume some order. The order is that, the value at current time is
	 * in first d (d = dimension() of the curve) rows of C, the the first jet past the current time is
	 * in the following d rows,..., finally last d rows is the Jet at pastTime().
	 * see also DoubletonInterface.
	 */
	BaseClass& set_Cr0(MatrixType const &C, VectorType const &r0);
	/** see DoubletonInterface */
	BaseClass& set_B(MatrixType const &B);
	/** see DoubletonInterface */
	BaseClass& set_r(VectorType const &r);
	/** see DoubletonInterface */
	BaseClass& set_Binv(MatrixType const &Binv);
	/** see DoubletonInterface */
	BaseClass& set_x(VectorType* x, bool passOwnership = false);
	/** see DoubletonInterface */
	BaseClass& set_C(MatrixType* C, bool passOwnership = false);
	/** see DoubletonInterface */
	BaseClass& set_r0(VectorType* r0, bool passOwnership = false);
	/** see DoubletonInterface */
	BaseClass& set_Cr0(MatrixType* C, VectorType* r0, bool passCOwnership = false, bool passr0Ownership = false);
	/** see DoubletonInterface */
	BaseClass& set_B(MatrixType* B, bool passOwnership = false);
	/** see DoubletonInterface */
	BaseClass& set_r(VectorType* r, bool passOwnership = false);
	/** see DoubletonInterface */
	BaseClass& set_Binv(MatrixType* B, bool passOwnership = false);
	/** see DoubletonInterface */
	VectorType* take_x();
	/** see DoubletonInterface */
	MatrixType* take_C();
	/** see DoubletonInterface */
	VectorType* take_r0();
	/** see DoubletonInterface */
	MatrixType* take_B();
	/** see DoubletonInterface */
	VectorType* take_r();
	/** see DoubletonInterface */
	MatrixType* take_Binv();
	using BaseClass::makeStorage_x;		///< see DoubletonInterface
	using BaseClass::makeStorage_C;		///< see DoubletonInterface
	using BaseClass::makeStorage_r0;	///< see DoubletonInterface
	using BaseClass::makeStorage_B;		///< see DoubletonInterface
	using BaseClass::makeStorage_r;		///< see DoubletonInterface
	using BaseClass::hull;  			///< see DoubletonInterface
	using BaseClass::midPoint;  		///< see DoubletonInterface
	/** see DoubletonInterface */
	BaseClass& add(VectorType const & v);
	/** see DoubletonInterface */
	BaseClass& add(BaseClass const & set);
	/** see DoubletonInterface */
	BaseClass& mulThenAdd(ScalarType const& c, BaseClass const & set);
	/** see DoubletonInterface */
	BaseClass& mul(ScalarType const& c);

	/** allows to set Xi for all Jets */
	BaseClass& set_Xi(VectorType const &Xi){
		size_type d = dimension();
		if (Xi.dimension() != d * m_pieces.size())
			throw std::logic_error("DDESolutionCUrve::set_Xi: incompatible dimension.");
		size_type I = 0;
		for (auto jet = rbegin(); jet != rend(); ++jet){
			VectorType Xi_part(d);
			for (size_type i = 0; i < d; ++i, ++I) Xi_part[i] = Xi[I];
			(*jet)->set_Xi(Xi_part);
		}
		return *this;
	}

	/** allows to retrieve Xi for all Jets */
	VectorType get_Xi() const{
		size_type d = dimension();
		VectorType Xi(d * m_pieces.size());
		size_type I = 0;
		for (auto jet = rbegin(); jet != rend(); ++jet){
			VectorType Xi_part = (*jet)->get_Xi();
			for (size_type i = 0; i < d; ++i, ++I) Xi[I] = Xi_part[i];
		}
		return Xi;
	}

	/** returns index in the array of the Jets for a given TimePoint */
	size_type pointToIndex(TimePointType const& t) const {
		int ibase = int(pastTime());
		int i0 = int(currentTime()) - ibase;
		int i = int(t) - ibase;
		if (i < 0 || i >= i0)
			throw std::domain_error("DDESolutionCurve::pointToIndex(): piece index outside range.");
		return i;
	}
	CurvePieceType& getPiece(TimePointType const& t) const {
		size_type i;
		try { i = pointToIndex(t); } catch (std::domain_error &e) { throw rethrow("DDESolutionCurve::getPiece(TimePoint):", e); }
		return *(m_pieces[i]);
	}
	/** it assumes that index is already processed (see pointToIndex()) */
	CurvePieceType& getPiece(size_type const& index) const {
		int ibase = int(pastTime());
		int i0 = int(currentTime()) - ibase;
		if (index < 0 || index >= i0)
			throw std::domain_error("DDESolutionCurve::getPiece(size_type): piece index outside range.");
		return *(m_pieces[index]);
	}
	VectorType eval(TimePointType t) const {
		try{
			return t == t0() ? VectorType(m_valueAtCurrent) : getPiece(t).eval(t);
		} catch (std::domain_error &e){
			throw rethrow("DDESolutionCurve::eval(TimePoint):", e);
		}
	}
	/**
	 * WARNING: THIS FUNCTION IS NOT RIGOROUS FOR INTERVALS! DO NOT USE IN PROOFS
	 * WARNING: COMPUTE RIGOROUSLY ONLY USING eval(TimePointType t) VERSION
	 */
	VectorType eval(RealType t) const {
		// TODO: (FUTURE???) raise warning (compilation time?) when this function is used?
		// TODO: it seems now the procudre looks safe, check!
		RealType epsi;
		TimePointType ti = m_grid.point(0);
		m_grid.split(t, ti, epsi);
		// std::cerr << "eval t = " << ti << " " << epsi << std::endl;
		try{
			return (ti == t0() && epsi == 0.0) ? VectorType(m_valueAtCurrent) : /* getPiece(ti).eval(t); */ getPiece(ti).evalAtDelta(epsi);
		} catch (std::domain_error &e){
			throw rethrow("DDESolutionCurve::eval(Real):", e);
		}
	}
	/** this is rigorous, it returns the value stored in the Jet at a given time point of the grid. */
	void eval(TimePointType t, SetType& out) const {
		try{
			out = (t == t0() ? m_valueAtCurrent : getPiece(t)[0]);
		} catch (std::domain_error &e){
			throw rethrow("DDESolutionCurve::eval(TimePoint, SetType&):", e);
		}
	}

	void eval(RealType t, SetType& out) const {
		throw std::logic_error("DDESolutionCurve::eval(RealType, SetType&): Not implemented yet.");
	}

	/*
	 * add piece of curve by pointer,
	 * optionally pass ownership of the pointer
	 * This is for optimization purposes.
	 * returns this for a cascade
	 * THIS IS FOR OPTIMIZATION PURPOSES AND SHOULD BE USED BY ADVENCED USER ONLY
	 */
	DDESolutionCurve& addPiece(CurvePieceType* newPiece, bool passOwnership = false);
	/*
	 * add piece of curve by making a copy,
	 * the copy will be owned by the DDESolutionCurve
	 * returns this for a cascade
	 * THIS SHOULD BE PREFFERED WAY TO ADD PIECES BY END USER
	 */
	DDESolutionCurve& addPiece(CurvePieceType const & newPiece);
	/** same as addPiece but adds in the past (before existing pieces) */
	DDESolutionCurve& addPastPiece(CurvePieceType* newPiece, bool passOwnership = false);
	/** same as addPiece but adds in the past (before existing pieces) */
	DDESolutionCurve& addPastPiece(const CurvePieceType& newPiece);
	/** set the value at a current time */
	DDESolutionCurve& setValueAtCurrent(const SetType& value);
	/** set the value at a current time */
	DDESolutionCurve& setValueAtCurrent(const VectorType& value);
	/** get the value at current */
	SetType& getValueAtCurrent() { return m_valueAtCurrent; };
	/** get the value at current (const object) */
	const SetType& getValueAtCurrent() const { return m_valueAtCurrent; };

	JetType jet(RealType t0, RealType t1) const { throw std::logic_error("DDESolutionCurve::jet: only supports jets over grid intervals."); }
	JetType jet(TimePointType t0, TimePointType t1) const {
		if (t0 + 1 != t1) throw std::logic_error("DDESolutionCurve::jet: does not support intervals longer than grid step size");
		try { return getPiece(t0); } catch (std::domain_error &e){ throw rethrow("DDESolutionCurve::jet(TimePoints):", e); }
	}
	JetType jet(TimePointType t0) const {
		try { return getPiece(t0); } catch (std::domain_error &e){ throw rethrow("DDESolutionCurve::jet(TimePoint):", e); }
	}
	/** use with caoution! Mainly for internal use for speed purposes */
	JetType* jetPtr(TimePointType t0) const {
		try { return &getPiece(t0); } catch (std::domain_error &e){ throw rethrow("DDESolutionCurve::jetPtr(TimePoint):", e); }
	}
	SetType j(size_type i, size_type k) const { try { return getPiece(i)[k]; } catch (std::domain_error &e){ throw rethrow("DDESolutionCurve::j(index, order): ", e); } }
	SetType j(TimePointType t0, size_type k) const { try { return getPiece(t0)[k]; } catch (std::domain_error &e){ throw rethrow("DDESolutionCurve::j(TimePoint, order): ", e); } }
	JetType const & j(TimePointType t0) const { try { return getPiece(t0); } catch (std::domain_error &e){ throw rethrow("DDESolutionCurve::j(TimePoint): ", e); } }
	JetType& j(TimePointType t0) { try { return getPiece(t0); } catch (std::domain_error &e){ throw rethrow("DDESolutionCurve::j(TimePoint): ", e); } }
	//SetType j(AnnotationType& annotation) const { try { return getPiece(annotation.i)[annotation.k]; } catch (std::domain_error &e){ throw rethrow("DDESolutionCurve::j(annotation): ", e); } }

	/** returns order of jet available at point t. This is for optimization purposes (do not copy whole Jet just to know its order) */
	size_type jetOrderAt(TimePointType t0){
		try { return getPiece(t0).order(); } catch (std::domain_error &e){ throw rethrow("DDESolutionCurve::jetorderAt(TimePoint): ", e); }
	}

	const GridType& getGrid() const { return m_grid; }
	const TimePointType getStep() const { return m_grid.point(1); }

	Class midCurve() const{
		Class result = *this; // TODO: (???) this is not so optimal to make a copy then once more rewrite items in this copy... REFACTOR IT!
		// make midCurve of all the CurvePieces
		for(auto i = result.begin(); i != result.end(); ++i)
			(**i) = (**i).midCurve();
		// make sure that the r0 vector is compatible with what CurvePieces demands after midCurve()
		result.deallocateR0();
		if (result.begin() != result.end())
			result.m_r0 = new VectorType((**result.begin()).storageN0());
		else
			result.m_r0 = new VectorType();  // 0-dim vector
		// then we update value at current time t0
		result.setValueAtCurrent(SetType(result.m_valueAtCurrent.midPoint(), result.m_r0));
		// and then update all r0 to be the same.
		result.updateCommonR0();
		return result;
	}
	/** compute derivative of the whole curve w.r.t. t of order k */
	Class dt(int k = 1) const{
		Class result = *this; // TODO: (???) this is not so optimal to make a copy then once more rewrite items in this copy... REFACTOR IT!
		for(auto i = result.begin(); i != result.end(); ++i)
			(**i) = (**i).dt(k);
		// TODO: IMPORTANT! what about the head????? do it as in Nonrig case? We should supply solver to move by one.
		return result;
	}

	/** show detailed info */
	std::string show() const;
	/** friend operator must be inline, or the linker gets confused */
	friend std::ostream& operator<<(std::ostream& out, DDESolutionCurve<SetSpec> const &  curve){ curve.writeTo(out); return out; }

	/**
	 * TODO: (FUTURE) this should be not templated by DynSysSpec but concrete "interfaced" type, known in advance
	 */
	template<typename DynSysSpec>
	DDESolutionCurve& move(DynSysSpec& solver);

	/**
	 * Naming convention. To be a pair of extend <-> epsilonShift
	 * See move() for more docs.
	 */
	template<typename DynSysSpec>
	DDESolutionCurve& extend(DynSysSpec& solver) { this->move(solver); }

	/**
	 * returns last enclosure of the value of the function over [t0, t0+h]
	 * obtained in extend() / move() procedure.
	 *
	 * TODO: (NOT URGENT) make a function that returns enclosure of the whole Jet?
	 */
	VectorType getLastEnclosure() const { return m_lastEnclosure; }

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
	 * The function does not check if the result is theoretically rigorous,
	 * i.e. it does not check if the solution curve is smooth enough to
	 * produce rigorous enclosure of a given order. Other classes
	 * of the library takes this responsibility (i.e. PoincareMap).
	 *
	 * TODO: (FUTURE) this should be not templated by DynSysSpec but concrete "interfaced" type, known in advance (?)
	 */
	template<typename DynSysSpec>
	void epsilonShift(
				DynSysSpec const& solver,
				RealType const& epsilon,
				DDESolutionCurve& out_result) const {
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
	 * The function does not check if the result is theoretically rigorous,
	 * i.e. it does not check if the solution curve is smooth enough to
	 * produce rigorous enclosure of a given order. Other classes
	 * of the library takes this responsibility (i.e. PoincareMap).
	 *
	 * TODO: (FUTURE) this should be not templated by DynSysSpec but concrete "interfaced" type, known in advance
	 */
	template<typename DynSysSpec>
	void epsilonShift(
				DynSysSpec const& solver,
				TimePointType const& at_t0, RealType const& epsilon,
				DDESolutionCurve& out_result) const;

	/** see base class DoubletonInterface */
	virtual void reinitialize(size_type d, size_type N0){
		throw std::logic_error("DDESolutionCurve::reinitialize: Not Implemented Yet");
	}

	/** computes dot taking appropriate number of first components from Curve */
	ScalarType dot(VectorType const& v) const;
	/** computes dot taking appropriate jets into account */
	template<typename JetSection>
	ScalarType dot(JetSection const& v) const;
//	/** this is L2 functional norm dot product of the two functions */
//	ScalarType dot_L2(DDESolutionCurve const& other) const;

	virtual ~DDESolutionCurve() { this->deallocate(); }

	/**
	 * it reduces the representation so that the maximal order
	 * of the jets is the same as in out_result.
	 *
	 * It throws exception if the out_result order is
	 * higher than this solution at some jet.
	 */
	void reduce(DDESolutionCurve& out_result) const;

protected:
	const GridType& m_grid;
	TimePointType m_t_current;
	VectorType m_lastEnclosure;
	/**
	 * this r0 will always be common for all elements.
	 * THIS IS IMPORTANT (because how the Lohner algorithm is implemented).
	 * Solution is always the owner of r0, all copies must copy it
	 * by creating new instance. Each instance must assure that
	 * m_r0 is always copied as a pointer (not owned) to the
	 * pieces and the current value. This variable must be before
	 * other dependent variables, because I use this fact
	 * to allocate it before others, so I can sometimes
	 * pass it to later dependable objects in the constructor list!
	 */
	VectorType* m_r0;
	SetType m_valueAtCurrent;
	size_type m_dimension;
	PiecesStorageType m_pieces;
	OwnerStorageType m_pieceOwner;

	// TODO: (NOT URGENT) clean-up the names and this whole allocation (e.g. allocateR0() etc. to be consistent)
	void allocate(size_type N0 = 0){ m_r0 = new VectorType(N0); }
	void deallocatePieces(){
		auto pIt = begin();
		auto oIt = m_pieceOwner.begin();
		for (; pIt != end(); ++pIt, ++oIt){
			if (*oIt) delete *pIt;
			// delete *pIt is ok here, as *pIt is a CurvePieceTypePtr
			// (**pIt is the CurvePiece object, sic!).
			*pIt = NULL;
		}
		m_pieces.clear();
		m_pieceOwner.clear();
	}
	void deallocateR0(){ capd::ddes::helper_safe_delete(m_r0, true); /* TODO: (RETHINK) r0 in pieces and valueAtCurrent is dangling now, should update to NULL i Think */ }
	void deallocate(){ deallocatePieces(); deallocateR0(); }
	void reallocateR0(size_type N0){ deallocateR0(); allocate(N0); }
	void reallocateR0(){ size_type N0 = storageN0(); deallocateR0(); allocate(N0); }
	void reallocate(size_type N0){ deallocate(); allocate(N0); }
	void reallocate(){ size_type N0 = storageN0(); deallocate(); allocate(N0); }
	void updateCommonR0(){
		try { m_valueAtCurrent.set_r0(m_r0); } catch (std::logic_error &e) {
			throw rethrow("DDESolutionCurve::updateCommonR0: at current value", e);
		}
		for (auto it = begin() ; it != end(); ++it){
			try {
				(*it)->set_r0(m_r0);
			} catch (std::logic_error &e) {
				throw rethrow("DDESolutionCurve::updateCommonR0: at some piece", e);
			}
		}
	}

	/** low level, does not care about changing time points, does not clear m_pieces, etc. just copies from to */
	template<typename IteratorSpec>
	void copyPieces(IteratorSpec from, IteratorSpec to){
		// TODO: (FUTURE) that's not the most optimal implementation...
		for (; from != to; ++from){
			m_pieces.push_back(new CurvePieceType(*from));
			m_pieceOwner.push_back(true);
		}
	}
	void copyPieces(iterator from, iterator to){
		// TODO: (NOT URGENT) that's not the most optimal implementation...
		for (; from != to; ++from){
			m_pieces.push_back(new CurvePieceType(**from));
			m_pieceOwner.push_back(true);
		}
	}
	void copyPieces(const_iterator from, const_iterator to){
		// TODO: (NOT URGENT) that's not the most optimal implementation...
		for (; from != to; ++from){
			m_pieces.push_back(new CurvePieceType(**from));
			m_pieceOwner.push_back(true);
		}
	}

	/**
	 * internal function to make short impl. of << operator here, and the long impl. in .hpp file
	 * (see comment in << operator)
	 * TODO: (NOT URGENT) maybe rethink... but I like how it is now with inline friend operator<<().
	 */
	void writeTo(std::ostream& out) const;
};


} // namespace ddes
} // namespace capd


#endif /* _CAPD_DDES_DDESOLUTIONCURVE_H_ */
