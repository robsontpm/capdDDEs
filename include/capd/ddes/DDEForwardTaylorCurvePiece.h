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

#ifndef _CAPD_DDES_DDEFORWARDTAYLORCURVEPIECE_H_
#define _CAPD_DDES_DDEFORWARDTAYLORCURVEPIECE_H_

#include <capd/ddes/DDECommon.h>
#include <capd/basicalg/TypeTraits.h>
#include <stdexcept>
#include <string>
#include <capd/ddes/storage/storage.h>
#include <capd/ddes/storage/GenericJet.h>

namespace capd{
namespace ddes{

/**
 * A class to represent a part of some curve $x$ over $[t0, t1)$ given by:
 *
 *    $(J^{[N]}_{t0}x)(t) + (n+1) * \int_{t0}^{t} \xi(s) (t0-s)^{N} ds$
 *
 * where $\xi(s) \in R \subset \mathbb{R}^d$
 *
 * in the code below:
 * d = dimension
 * n = order
 *
 * If the interval version is used, then the class can represent
 * sets of Curves, so that $J^{[N]}_{t0} \in x + A * r0 + Q * r$,
 * (x is here a middle point of the set and Jet is interpreted as
 * a collection of coefficients).
 * $r0$ is some vector of any given dimension, compatible with $A$.
 * We denote the dimension of $r0$ as $N0$.
 * $Q$ is a diagonal block matrix, so that $Q_ii$ is orthonormal $d \times d$ 
 * matrix and and $Q$ is 0 elsewhere. Vector $r$ is of dimension $d * n$.
 * we store matrices row-wise, ie in the code below
 * m_A[k] represents a $d \times N0$  matrix and m_Q[k] is $Q_kk$ a 
 * $d \times d$ matrix. Finally, we store $r$ and $x$ again component-wise, 
 * so that in code below m_r[k] and $m_x[k]$ is a $d$-dimensional vector.
 * Therefore we have: $(J^{[N]}_{t0}x)_{[k]} \in m_x[k] + m_A[k] * r0 + m_Q[k] * m_r[k]$
 *
 * NOTE: we store $m_r0$ as a one $N0$-dimensional vector, but for optimization
 *       purposes we hold the variable as a pointer - to allow sharing of this 
 *       vector among many instances (important in DDESolutionCurve)
 *
 *
 * TODO: (NOT URGENT) extract abstract interface that has eval, taylor, summa, etc?
 * TODO: (NOT URGENT) base interface is already done as GenericJet, it has taylor (or eval).
 * TODO: (NOT URGENT) This should inherit Jet structure from there and add summa and taylor, and reload eval.
 */
template<
	typename TimePointSpec,
	typename SetSpec,
	bool isInterval = capd::TypeTraits<typename SetSpec::MatrixType::ScalarType>::isInterval
>
class DDEForwardTaylorCurvePiece :
		public DoubletonInterface<typename SetSpec::MatrixType, typename SetSpec::Policy>{
public:
	// TODO: (FAR FUTURE, NOT URGENT) MatrixType (and others) must be of dynamic size (do some TraitsType check?)
	typedef SetSpec DoubletonStorageType;
	typedef SetSpec SetType;
	typedef typename SetType::MatrixType MatrixType;
	typedef typename SetType::Policy Policy;
	typedef typename SetType::QRPolicy QRPolicy;
	typedef DDEForwardTaylorCurvePiece<TimePointSpec, SetSpec, isInterval> Class;
	typedef DoubletonInterface<MatrixType, Policy> BaseClass;
	typedef typename SetType::VectorType VectorType;
	typedef typename SetType::ScalarType ScalarType;
	typedef typename SetType::size_type size_type;
	typedef ScalarType RealType;  // TODO: (FAR FUTURE, NOT URGENT) eventually do TypeTraits like Daniel's code
	typedef TimePointSpec TimePointType;
	typedef Class JetType;

	typedef const DoubletonStorageType* const_iterator;
	typedef DoubletonStorageType* iterator;

	typedef capd::ddes::GenericJet<TimePointType, VectorType, VectorType, MatrixType> ExternalJetType;

	/** should be one word */
	static std::string badge() { return "DDEForwardTaylorCurvePiece"; }

	/** copy operator. Works as in case of copy constructor (see doc there) */
	DDEForwardTaylorCurvePiece& operator=(DDEForwardTaylorCurvePiece const & orig);
	/** creates an empty curve of order 0 and dimension 0, by default at the point $t_0=0$ */
	DDEForwardTaylorCurvePiece(TimePointType t0 = TimePointType());
	/**
	 * copy constructor. if orig has a shared r0 and/or Xi then it will share it with the copy (this)
	 * Otherwise, it assures that the r0 / Xi is copied accordingly and the jet updated if necessary.
	 */
	DDEForwardTaylorCurvePiece(DDEForwardTaylorCurvePiece const & orig);
	/** makes a constant function of the vector value 0 */
	DDEForwardTaylorCurvePiece(TimePointType t0, size_type dimension, size_type order);
	/** makes a constant function of the vector value v */
	DDEForwardTaylorCurvePiece(TimePointType t0, size_type order, const VectorType& v, size_type N0 = 0);
	/** makes a constant function of the vector value v */
	DDEForwardTaylorCurvePiece(TimePointType t0, size_type order, const VectorType& v, VectorType* r0, bool passOwnership = false);
	/** fills vector with data. Determine order from the size of data. */
	DDEForwardTaylorCurvePiece(TimePointType t0, VectorType const* it, VectorType const* itEnd, size_type N0 = 0);
	/** fills vector with data. Determine order from the size of data. */
	DDEForwardTaylorCurvePiece(TimePointType t0, std::vector<VectorType> coeffs, size_type N0 = 0);
	/** makes a curve piece of constant function, located at t0, of a given order and of given value. The r0 will be taken from the sets r0 (and no ownership transfer occurs). */
	DDEForwardTaylorCurvePiece(TimePointType t0, size_type order, SetType& value);
	/** makes a curve piece of constant function, located at t0, of a given order and of given value. The r0 will be taken from the sets r0 (as copy). */
	DDEForwardTaylorCurvePiece(TimePointType t0, size_type order, const SetType& value);
	/** almost complete constructor, setups everything from scratch (computes invB) */
	DDEForwardTaylorCurvePiece(
			TimePointType t0,
			std::vector<VectorType> coeffs,
			std::vector<MatrixType> Cs,
			VectorType* r0, bool passR0Ownership,
			std::vector<MatrixType> Bs,
			std::vector<VectorType> rs,
			VectorType* Xi, bool passXiOwnership);
	/** most complete constructor, setups everything from scratch */
	DDEForwardTaylorCurvePiece(
			TimePointType t0,
			std::vector<VectorType> coeffs,
			std::vector<MatrixType> Cs,
			VectorType* r0, bool passR0Ownership,
			std::vector<MatrixType> Bs,
			std::vector<MatrixType> invBs,
			std::vector<VectorType> rs,
			VectorType* Xi, bool passXiOwnership);

	/** dimension of the value space of the Jet (i.e. $J : \R \to \R^{dimension}$) */
	size_type dimension() const { return m_dimension; }
	/** order of the Taylor part of the Jet */
	size_type order() const { return m_order; }

	/** DoubletonStorageInterface: returns a dimension of the Jet as a Vector (sequence) to store all the coefficients (without Xi part) */
	size_type storageDimension() const { return m_dimension * (m_order + 1); }
	/** DoubletonStorageInterface: returns a dimension of the r0 vector */
	size_type storageN0() const { return m_r0->dimension(); }
	/** see doubleton interface */
	VectorType get_x() const;
	/** see doubleton interface */
	MatrixType get_C() const;
	/** see doubleton interface */
	MatrixType get_B() const;
	/** see doubleton interface */
	VectorType get_r() const;
	/** see doubleton interface */
	VectorType get_r0() const;
	/** see doubleton interface */
	BaseClass& set_x(VectorType const &x);
	/** see doubleton interface */
	BaseClass& set_C(MatrixType const &C);
	/** see doubleton interface */
	BaseClass& set_r0(VectorType const &r0);
	/** see doubleton interface */
	BaseClass& set_Cr0(MatrixType const &C, VectorType const &r0);
	/**
	 * see doubleton interface
	 *
	 * WARNING: we assume additionaly, that B has a block-diagonal form, with d x d diagonal blocks (d = dimension()).
	 *          the code will raise exception if off-block-diagonal element is found non zero.
	 */
	BaseClass& set_Binv(MatrixType const &Binv);
	/**
	 * see doubleton interface
	 *
	 * WARNING: we assume additionaly, that B has a block-diagonal form, with d x d diagonal blocks (d = dimension()).
	 *          the code will raise exception if off-block-diagonal element is found non zero.
	 */
	BaseClass& set_B(MatrixType const &B);
	/** see doubleton interface */
	BaseClass& set_r(VectorType const &r);
	using BaseClass::makeStorage_x;		///< see Doubleton interface
	using BaseClass::makeStorage_C;		///< see Doubleton interface
	using BaseClass::makeStorage_r0;	///< see Doubleton interface
	using BaseClass::makeStorage_B;		///< see Doubleton interface
	using BaseClass::makeStorage_r;		///< see Doubleton interface
	using BaseClass::hull;				///< see Doubleton interface
	using BaseClass::midPoint; 			///< see Doubleton interface
	using BaseClass::dot;			 	///< see Doubleton interface

	/** multiply set by a scalar. Should take set structure into consideration. */
	virtual Class& mul(ScalarType const & c){
		for (auto jk = beginJet(); jk != endJet(); ++jk) jk->mul(c);
		(*this->m_Xi) *= c;
		return *this;
	}
	/** see doubleton interface */
	BaseClass& affineTransform(MatrixType const &M, VectorType const &v) { throw std::logic_error("DDEForwardTaylorCurvePiece::affineTransform: Not implemented yet"); }
	/** see doubleton interface */
	BaseClass& translate(VectorType const &v) { throw std::logic_error("DDEForwardTaylorCurvePiece::translate: Not implemented yet"); }

	/**
	 * computes dot product of this with JetType (does not take Xi part into consideration!)
	 *
	 * This is a simple dot product in vector space. Consider implementing some other like integral
	 * (but it should be equivalent more or less becouse of polynomials)
	 */
	ScalarType dot(ExternalJetType const& s) const {
		ScalarType value;
		auto ithis = this->beginJet();
		for (auto is = s.begin(); is != s.end(); ++is, ++ithis){
			value += ithis->dot(*is);
		}
		return value;
	}

	/** returns Forward Taylor Jet of the derivative of order n w.r.t. time, n must be smaller than the order */
	DDEForwardTaylorCurvePiece dt(size_type n = 1) const;
	/** evaluates the jet at time t0 + delta_t , delta_t should be positive!  */
	VectorType evalAtDelta(const RealType& delta_t) const;
	/** evaluates the jet at time t (please note that it is true time, not some small h w.r.t. jet.t0().  */
	VectorType eval(const RealType& t) const { return this->evalAtDelta(t - RealType(t0())); };
	/** helper. Evaluates n-th derivative (more precisely j^[n] = j^(n)/n! - it is used in the algorithm most often) of this jet at t. It should be faster than jet.dt(n).evalAtDelta(t)/n!. The delta_t should be positive! */
	VectorType evalCoeffAtDelta(size_type n, const RealType& delta_t) const;
	/** helper. Evaluates n-th derivative (more precisely j^[n] = j^(n)/n! - it is used in the algorithm most often) of this jet at t. It should be faster than jet.dt(n).eval(t)/n! */
	VectorType evalCoeff(size_type n, const RealType& t) const { return this->evalCoeffAtDelta(n, t - RealType(t0())); }
	/** evaluates the Taylor part of the jet at t0 + delta_t */
	VectorType taylorAtDelta(const RealType& delta_t) const;
	/** evaluates the Taylor part of the jet at t */
	VectorType taylor(const RealType& t) const { return this->taylorAtDelta(t - RealType(t0())); }
	/** evaluates the integral of Xi part of the jet at t0 + delta_t */
	VectorType summaAtDelta(const RealType& delta_t) const;
	/** evaluates the integral of Xi part of the jet at t */
	VectorType summa(const RealType& t) const { return this->summaAtDelta(t - RealType(t0())); }
	/** returns jet at t for a function represented by this. It is rigorous, so it returns validated estimates (overestimates) */
	JetType jetAt(const TimePointType& t) const { if (t == m_t0) { return *this; } else { return this->jetAt((RealType)t); } }
	/** returns jet at t for a function represented by this. It is rigorous, so it returns validated estimates (overestimates) */
	JetType jetAt(const RealType& t) const { throw std::logic_error("DDEForwardTaylorCurvePiece::jetAt(RealType t): Not implemented yet"); };
	/**
	 * evaluates the Taylor part of the jet at t0 + delta_t, delta_t should be positive.
	 * Note: out should be set representing [0,0]^d. Procedure do not test for it!
	 */
	void taylorAtDelta(const RealType& delta_t, SetType& out) const;
	/**
	 * evaluates the Taylor part of the jet at t, delta_t should be positive.
	 * * Note: out should be set representing [0,0]^d. Procedure do not test for it!
	 */
	void taylor(const RealType& t, SetType& out) const { this->taylorAtDelta(t - RealType(t0()), out); }
	/**
	 * evaluates the jet at time t0 + delta_t, delta_t should be positive.
	 * Note: out should be set representing [0,0]^d. Procedure do not test for it!
	 */
	void evalAtDelta(const RealType& delta_t, SetType& out) const;
	/**
	 * evaluates the jet at time t0 + delta_t, delta_t should be positive!
	 * Note: out should be set representing [0,0]^d. Procedure do not test for it!
	 */
	void eval(const RealType& t, SetType& out) const { evalAtDelta(t - RealType(t0()),  out); } ;
	/**
	 * evaluates the jet at time t0 + delta_t, delta_t should be positive!
	 * * Note: out should be set representing [0,0]^d. Procedure do not test for it!
	 */
	void evalCoeffAtDelta(size_type n, const RealType& delta_t, SetType& out) const;
	/**
	 * evaluates the jet at time t0 + delta_t, delta_t should be positive!
	 * Note: out should be set representing [0,0]^d. Procedure do not test for it!
	 */
	void evalCoeff(size_type n, const RealType& t, SetType& out) const { evalDtAtDelta(n, t - RealType(t0()),  out); } ;

	/** returns jet with Taylor part of diameter 0 and the same Xi */
	DDEForwardTaylorCurvePiece midCurve() const;
	/** returns true, if the diamater of Taylor part is 0 */
	bool isMidCurve() const { return storageN0(); }
	/** more verbose output, more human-friendly. Not rigorous-friendly. */
	std::string show() const;

	/** iterator to the 0-th order Taylor coefficient (might by many-dimensions) */
	iterator beginJet() { return m_jet_at_t0; }
	/** iterator to the past the order place. Standard thing in C++ */
	iterator endJet() { return m_jet_at_t0 + m_order + 1; }
	/** iterator to the last element in coefficients. Standard thing in C++ */
	iterator backJet() { return m_jet_at_t0 + m_order; }
	/** iterator to the 0-th order Taylor coefficient (might by many-dimensions) */
	const_iterator beginJet() const { return m_jet_at_t0; }
	/** iterator to the past the order place. Standard thing in C++ */
	const_iterator endJet() const { return m_jet_at_t0 + m_order + 1; }
	/** iterator to the last element in coefficients. Standard thing in C++ */
	const_iterator backJet() const { return m_jet_at_t0 + m_order; }
	/** returns n-th order Taylor coefficient. */
	DoubletonStorageType& operator[](size_type k){ return m_jet_at_t0[k]; }
	/** returns n-th order Taylor coefficient. */
	DoubletonStorageType const& operator[](size_type k) const { return m_jet_at_t0[k]; }

	/** set this vector to represent a constant function f(t) = value for all t > t0, the r0 vector is taken from the value set. */
	Class& setAsConstant(SetType& value){
		deallocateR0(); // this makes some pointers dangling (if set earlier), we will correct it later (when setting the value for all Jets).
		this->m_r0 = value.take_r0();
		this->m_r0_owner = true;
		VectorType zero = VectorType(value) * ScalarType(0);
		m_jet_at_t0[0] = value;
		for (auto jetIt = beginJet() + 1; jetIt != endJet(); ++jetIt)
			*jetIt = DoubletonStorageType(zero, this->m_r0);
		deallocateXi();
		set_Xi(zero);
		// no need to updateCommonR0(); - we have done it before (in set constructors)
		return *this;
	}
	/** set this vector to represent a constant function f(t) = value for all t > t0, the r0 vector is taken from the value set. */
	Class& setAsConstant(const SetType& value, VectorType* overwrite_r0 = NULL, bool passOwnership = false){
		deallocateR0(); // this makes some pointers dangling, we will correct it later (when setting the value for all Jets).
		this->m_r0 = overwrite_r0 ? overwrite_r0 : new VectorType(value.get_r0());
		this->m_r0_owner = overwrite_r0 ? passOwnership : true;
		VectorType zero = VectorType(value) * ScalarType(0);
		m_jet_at_t0[0] = value;
		for (auto jetIt = beginJet() + 1; jetIt != endJet(); ++jetIt)
			*jetIt = DoubletonStorageType(zero, this->m_r0);
		set_Xi(zero);
		// no need to updateCommonR0(); - we have done it before (in set constructors)
		return *this;
	}
	/** set this vector to represent a constant function f(t) = value for all t > t0 */
	Class& setAsConstant(VectorType const & value, size_type N0 = 0){
		reallocateR0(N0); // sets R0 to N0-dim zero vector.
		VectorType zero = value * ScalarType(0);
		m_jet_at_t0[0] = DoubletonStorageType(value, this->m_r0);
		for (auto jetIt = beginJet() + 1; jetIt != endJet(); ++jetIt)
			*jetIt = DoubletonStorageType(zero, this->m_r0);
		set_Xi(zero);
		// no need to updateCommonR0(); - we have done it before (in set constructors)
		return *this;
	}
	/** set this vector to represent a constant function f(t) = value for all t > t0 */
	Class& setAsConstant(VectorType const & value, VectorType* r0, bool passOwnership = false){
		deallocateR0();
		this->m_r0 = r0;
		this->m_r0_owner = passOwnership;
		VectorType zero = value * ScalarType(0);
		m_jet_at_t0[0] = DoubletonStorageType(value, this->m_r0);
		for (auto jetIt = beginJet() + 1; jetIt != endJet(); ++jetIt)
			*jetIt = DoubletonStorageType(zero, this->m_r0);
		set_Xi(zero);
		// no need to updateCommonR0(); - we have done it before (in set constructors)
		return *this;
	}
	/** set the base time of the jet */
	Class& setT0(TimePointType const & t0){ m_t0 = t0; return *this; }
	/** returns time at which this jet is located */
	TimePointType getT0() const { return m_t0; }
	/** naming convention. See getT0() */
	TimePointType t0() const { return m_t0; }

	/** set the estimate on the \xi part of the jet, makes object to own it. */
	void set_Xi(VectorType const& xi){ dimCheck(&xi); reallocateXi(); *m_Xi = xi; m_Xi_owner = true; }
	/** set the estimate on the \xi part of the jet, by default id does not pass ownership. Throws exception if Xi is of bad dimension. */
	void set_Xi(VectorType* xi, bool passOwnership = false){ dimCheck(xi); deallocateXi(); m_Xi = xi; m_Xi_owner = passOwnership; }
	/** returns a value of Xi */
	VectorType get_Xi() const { return *m_Xi; }
	/** returns Xi and forgets that object is owner (user is responsible for deallocating Xi) */
	VectorType* take_Xi() { this->m_Xi_owner = false; return m_Xi; }
	/** returns r0 and forgets that object is owner (user is responsible for deallocating r0) */
	VectorType* take_r0() { this->m_r0_owner = false; return m_r0; }
	/** set the estimate on the r0 part of the set, by default id does not pass ownership. Updates jet pointers! Throws exception if r0 is of bad dimension. */
	Class& set_r0(VectorType* r0, bool passOwnership = false) {
		if (!r0)
			throw std::logic_error("DDEForwardTaylorCurvePiece::set_r0(ptr): ptr is null!");
		if (r0->dimension() != storageN0()){
			std::ostringstream info;
			info << "DDEForwardTaylorCurvePiece::set_r0(ptr): vector has incompatible N0 dimension, ";
			info << "is " << r0->dimension() << ", expected " << storageN0() << ".";
			throw std::logic_error(info.str());
		}
		deallocateR0(); m_r0 = r0;
		m_r0_owner = passOwnership;
		updateCommonR0();
		return *this;
	}

	/** standard thing */
	virtual ~DDEForwardTaylorCurvePiece() { this->deallocate(); };

	/** see base class */
	virtual void reinitialize(size_type d, size_type N0){
		throw std::logic_error("DDEForwardTaylorCurvePiece::reinitialize: Not Supported Yet.");
	}

	friend std::ostream& operator<<(std::ostream & out, DDEForwardTaylorCurvePiece const & jet) {
		out << DDEForwardTaylorCurvePiece::badge() << std::endl;
		out << jet.m_dimension << " " << jet.m_order << std::endl;
		out << jet.getT0() << std::endl;
		out << "CommonR0" << std::endl;
		out << jet.get_r0() << std::endl;
		out << "Jet" << std::endl;
		size_type k = 0;
		for (auto kth = jet.beginJet(); kth != jet.endJet(); ++kth, ++k){
			out << k << "-th" << std::endl;
			out << (*kth).get_x() << std::endl;
			out << (*kth).get_C() << std::endl;
			out << (*kth).get_B() << std::endl;
			out << (*kth).get_Binv() << std::endl;
			out << (*kth).get_r() << std::endl;
		}
		out << "Xi" << std::endl;
		out << *(jet.m_Xi) << std::endl;
		return out;
	}

	friend std::istream& operator>>(std::istream & in, DDEForwardTaylorCurvePiece & jet) {
		helper_dump_badge(in);
		jet.deallocate();
		in >> jet.m_dimension >> jet.m_order;
		in >> jet.m_t0;
		helper_dump_badge(in); // "CommonR0"
		VectorType r0;
		in >> r0;
		jet.allocate(r0.dimension());
		(*jet.m_r0) = r0;
		MatrixType C, B, Binv;
		VectorType x, r;
		helper_dump_badge(in); // "Jet"
		for (size_type k = 0; k <= jet.m_order; ++k){
			helper_dump_badge(in); // "k-th"
			in >> x;
			in >> C;
			in >> B;
			in >> Binv;
			in >> r;
			jet.m_jet_at_t0[k] = DoubletonStorageType(x, C, jet.m_r0, B, Binv, r);
		}
		helper_dump_badge(in); // "Xi"
		in >> (*jet.m_Xi);
		return in;
	}

protected:
	TimePointType m_t0;
	size_type m_dimension;
	size_type m_order;
	DoubletonStorageType* m_jet_at_t0;
	VectorType* m_r0;
	bool m_r0_owner;
	VectorType* m_Xi;
	bool m_Xi_owner;

	/** checks if the dimension of v is compatible with current set */
	void dimCheck(const VectorType* const v){
		if (!v || v->dimension() != dimension()){
			std::ostringstream info;
			info << "DDEForwardTaylorCurvePiece::dimCheck(): vector has incompatible dimension: ";
			info << "was " << v->dimension() << ", expected " << m_dimension << ".";
			throw std::logic_error(info.str());
		}
	}

	/** checks if the dimension of C matrix is compatible with current set C matrix (might be swapped) */
	void dimCheckC(const MatrixType* const C){
		if (!C || C->numberOfRows() != dimension() || C->numberOfColumns() != storageN0()){
			std::ostringstream info;
			info << "DDEForwardTaylorCurvePiece::dimCheck(): C-part matrix has incompatible dimension: ";
			info << "was M(" << C->numberOfRows() << "," << C->numberOfColumns();
			info << ", expected M(" << dimension() << "," << storageN0() << ".";
			throw std::logic_error(info.str());
		}
	}
	/** checks if the dimension of B is compatible with current set B matrix */
	void dimCheckB(const MatrixType* const B){
		if (!B || B->numberOfRows() != dimension() || B->numberOfColumns() != dimension()){
			std::ostringstream info;
			info << "DDEForwardTaylorCurvePiece::dimCheck(): B-part matrix has incompatible dimension: ";
			info << "was M(" << B->numberOfRows() << "," << B->numberOfColumns();
			info << ", expected M(" << dimension() << "," << dimension() << ".";
			throw std::logic_error(info.str());
		}
	}

	/** technical, manages deallocation of resources */
	void deallocateJet(){ capd::ddes::helper_safe_array_delete(m_jet_at_t0, true); }
	/** technical, manages deallocation of resources */
	void deallocateXi(){ capd::ddes::helper_safe_delete(m_Xi, m_Xi_owner); }
	/** technical, manages deallocation of resources */
	void deallocateR0(){ capd::ddes::helper_safe_delete(m_r0, m_r0_owner); /* m_jet can have dangling pointers now */ }
	/** technical, manages deallocation of resources */
	void deallocate(){  deallocateJet(); deallocateXi(); deallocateR0(); }
	/** technical, manages deallocation of resources */
	void allocateJet(){
		m_jet_at_t0 = new DoubletonStorageType[order()+1];
		VectorType zero(dimension());
		for (auto k = beginJet(); k != endJet(); ++k){
			if (m_r0)
				*k = DoubletonStorageType(zero, m_r0);
			else
				*k = DoubletonStorageType(zero, 0);
		}
	}
	/** technical, manages allocation of resources */
	void allocateXi(){ m_Xi_owner = true; m_Xi = new VectorType(dimension()); }
	/** technical, manages allocation of resources */
	void allocateR0(size_type N0 = 0){ m_r0_owner = true; m_r0 = new VectorType(N0); if (m_jet_at_t0) updateCommonR0(); }
	/** technical, manages allocation of resources */
	void allocate(size_type N0 = 0){ allocateR0(N0); allocateJet(); allocateXi(); updateCommonR0(); }
	/** technical, manages allocation of resources */
	void reallocateJet(){ deallocateJet(); allocateJet(); }
	/** technical, manages re-allocation of resources */
	void reallocateXi(){ deallocateXi(); allocateXi(); }
	/** technical, manages re-allocation of resources */
	void reallocateR0(size_type N0 = 0){ deallocateR0(); allocateR0(N0); updateCommonR0(); }
	/** technical, manages re-allocation of resources */
	void reallocate(size_type N0 = 0){ deallocate(); allocate(N0); }
	/** technical, takes care that all jets share the same r0 */
	void updateCommonR0(){ for (auto jetIt = beginJet() ; jetIt != endJet(); ++jetIt) jetIt->set_r0(m_r0); }

	/** setup data from a generic form of Iterator (to collection of Vectors) */
	template<typename IteratorType>
	void setupFromData(IteratorType it, IteratorType itEnd){
		if (itEnd <= it) throw std::logic_error("DDEForwardTaylorCurvePiece::copyFromData(): cannot copy from empty range.");
		m_order = itEnd - it - 1;
		m_dimension = it->dimension();  // we assume dimension of the first element
		reallocate(storageN0());					// this will unnecessary reallocate R0, rethink...
		MatrixType eachC(m_dimension, storageN0()); // to make sets of good shape to N0. The matrix = 0, so it should be set later by end user.
		std::cout << "storageN0-piece" << storageN0() << std::endl;
		for (auto jk = beginJet(); it != itEnd; ++it, ++jk){
			try{ dimCheck(&(*it)); } catch (std::logic_error& e) { throw rethrow("DDEForwardTaylorCurvePiece::copyFromData(): incompatible dimension.", e); }
			try{*jk = SetType(*it, eachC, *m_r0);} catch (std::logic_error& e) { throw rethrow("DDEForwardTaylorCurvePiece::copyFromData(): SetType constructor error.", e); }
		}
		// previously, in SetType(..., eachC, *m_r0) we used a copy of the m_r0. Now we make it shared.
		try{ updateCommonR0(); } catch (std::logic_error& e) { throw rethrow("DDEForwardTaylorCurvePiece::copyFromData(): updateCommonR0 error.", e); }
	}

	/**
	 * Setup set from a rich set of data.
	 * Please note that you pass r0 and Xi as pointers, but not pass ownership flags
	 * You need to manage flags manually, outside the procedure.
	 * This is highly internal function, not meant to be used by the end-user.
	 */
	template<typename XItSpec, typename CItSpec, typename BItSpec, typename RItSpec>
	void setupFromData(
			XItSpec xit, XItSpec xend,
			CItSpec Cit, CItSpec Cend,
			VectorType* r0,
			BItSpec Bit, BItSpec Bend,
			BItSpec invBit, BItSpec invBend,
			RItSpec rit, RItSpec rend,
			VectorType* Xi)
	{
		deallocate();
		if (xend <= xit || Cend <= Cit || Bend <= Bit || rend <= rit)
			throw std::logic_error("DDEForwardTaylorCurvePiece::copyFromData(): cannot copy from empty range.");
		m_order = xend - xit - 1;
		if (Cend - Cit <= m_order || Bend - Bit <= m_order || rend - rit <= m_order)
			throw std::logic_error("DDEForwardTaylorCurvePiece::copyFromData(): not enough data to construct object.");
		m_dimension = xit->dimension();
		m_Xi = Xi;
		m_r0 = r0;
		// we only allocate Jet, as r0 and Xi are already present as pointers.
		allocateJet();
		for (auto jk = beginJet(); jk != endJet(); ++jk, ++xit, ++Cit, ++Bit, ++invBit, ++rit){
			try{ dimCheck(&(*xit)); } catch (std::logic_error& e) { throw rethrow("DDEForwardTaylorCurvePiece::copyFromData(): incompatible dimension in x.", e); }
			try{ dimCheck(&(*rit)); } catch (std::logic_error& e) { throw rethrow("DDEForwardTaylorCurvePiece::copyFromData(): incompatible dimension in r.", e); }
			try{ dimCheckC(&(*Cit)); } catch (std::logic_error& e) { throw rethrow("DDEForwardTaylorCurvePiece::copyFromData(): incompatible dimension in C.", e); }
			try{ dimCheckB(&(*Bit)); } catch (std::logic_error& e) { throw rethrow("DDEForwardTaylorCurvePiece::copyFromData(): incompatible dimension in B.", e); }
			try{ dimCheckB(&(*invBit)); } catch (std::logic_error& e) { throw rethrow("DDEForwardTaylorCurvePiece::copyFromData(): incompatible dimension in invB.", e); }
			*jk = SetType(*xit, *Cit, r0, *Bit, *invBit, *rit);
		}
		// no need for updateCommonR0 - creation of sets already done that
	}

};


} // namespace ddes
} // namespace capd


#endif /* _CAPD_DDES_DDEFORWARDTAYLORCURVEPIECE_H_ */
