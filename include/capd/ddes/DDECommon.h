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

#ifndef _CAPD_DDES_DDECOMMON_H_
#define _CAPD_DDES_DDECOMMON_H_

#include <capd/capdlib.h>

#include <vector>		// I use std::vectors for convient containers in some places
#include <string>		// I pass for now exceptions with build-in std::logic_error mostly
#include <stdexcept>	// I use strongs in many routines for this purpose

namespace capd{
namespace ddes{

/** this is used in many places, so DRY */
template<typename PtrT>
void helper_safe_delete(PtrT& item, bool is_owner){
	if (item && is_owner){
		delete item;
		item = 0;
	}
}

/** this is used in many places, so DRY */
template<typename PtrT>
void helper_safe_array_delete(PtrT& item, bool is_owner){
	if (item && is_owner){
		delete[] item;
		item = 0;
	}
}

/** After reconsider, this probably will be removed. TODO: (NOT URGENT) probably remove */
void helper_dump_line(std::istream & in);

/** After reconsider, this probably will be removed. TODO: (NOT URGENT) probably remove  */
void helper_dump_badge(std::istream & in);

/** helper function, to be used to simply override versions for double (numerical nonrigorous) and Interval (rigorous) */
template<typename RealSpec>
RealSpec ecloseStep(RealSpec h){ return RealSpec(0.0, 1.0) * h; }
/** nonrigorous version for double numerical only method */
double ecloseStep(double h);

/** helper function for pretty printing. TODO: (NOT URGENT) move to other files (helper?) */
template<typename RealSpec>
std::string showEnclosedInterval(RealSpec a, RealSpec b){ std::ostringstream oss; oss << "[" << a.leftBound() << ", " << b.rightBound() << ")"; return oss.str(); }
/** nonrigorous version for double numerical only method TODO: (NOT URGENT) move to other files (helper?) */
std::string showEnclosedInterval(double a, double b);

/**
 * makes an exception with extra message in front and glued with glue. Default glue = "\n    "
 * makes a good work to produce some kind of stack trace.
 * Sugested usage:
 * try {
 *    ...
 * } catch (some_exception& e){
 *    throw rethrow("Some extra info: ", e);
 * }
 * We could make this function to make throw, but then it confuses
 * IDEs in code inspection, especially, when used in functions supposed
 * to return value in all execution paths.
 */
template<typename ExceptionSpec>
ExceptionSpec rethrow(std::string msg, ExceptionSpec const &e, std::string glue = "\n    ") {
	return ExceptionSpec(msg + glue + e.what());
}

/**
 * makes an exception with extra message in front and glued with glue. Default glue = "\n    "
 * makes a good work to produce some kind of stack trace.
 * Sugested usage:
 * try {
 *    ...
 * } catch (some_exception& e){
 *    throw rethrow("Some extra info: ", e);
 * }
 * We could make this function to make throw, but then it confuses
 * IDEs in code inspection, especially, when used in functions supposed
 * to return value in all execution paths.
 */
template<typename ExceptionSpec, typename NewExceptionSpec>
NewExceptionSpec rethrow(std::string msg, ExceptionSpec const &e, std::string glue = "\n    ") {
	return NewExceptionSpec(msg + glue + e.what());
}

/** returns closest int */
template<typename AnythingSpec>
int closestInt(AnythingSpec const & value){ return int(value); }
/** special case for capd::intervals::Interval template */
template<typename T_Bound, typename T_Rnd>
int closestInt(capd::intervals::Interval<T_Bound, T_Rnd> const & value){ return closestInt(value.leftBound()); }
/** returns closest int smaller than the value */
template<typename AnythingSpec>
int closestSmallerInt(AnythingSpec const & value){ return int(value) - (value < 0 ? 1 : 0); }
/** special case for capd::intervals::Interval template */
template<typename T_Bound, typename T_Rnd>
int closestSmallerInt(capd::intervals::Interval<T_Bound, T_Rnd> const & value){ return closestSmallerInt(value.leftBound()); }
/** special case for capd::interval template */
int closestSmallerInt(capd::interval const & value);

/** require C++11 */
template<typename VectorType, typename ...VectorTypes>
typename VectorType::size_type
commonDimension(VectorType const & first, VectorTypes const & ...others) {
    return first.dimension(sizeof...(others) ? commonDimension(others...) : 0);
}
/** require C++11 TODO: (NOT URGENT) probably remove  */
template<typename OutIterator, typename VectorType, typename ...VectorTypes>
void joinVectorsRaw(OutIterator& out, VectorType const & first, VectorTypes const & ...others) {
	for (auto vIt = first.begin(); vIt != first.end(); ++vIt, ++out)
		*out = *vIt;
	if (sizeof...(others))
		joinVectorsRaw(out, others...);
}
/** require C++11 TODO: (NOT URGENT) probably remove  */
template<typename VectorType, typename ...VectorTypes>
VectorType joinVectors(VectorType const & first, VectorTypes const & ...others) {
	typename VectorType::size_type dim = commonDimension(first, others...);
	VectorType result(dim);
	joinVectorsRaw(result.begin(), first, others...);
	return result;
}

/**
 * represents time points on a grid i * h, i - integer, h - real
 *
 * TODO: (NOT URGENT, FAR FUTURE, RETHINK) rewrite so that TimePoint holds ref to Grid, then no problems will be presented!
 * TODO: (NOT URGENT, FAR FUTURE, RETHINK) (one problem: when we want to have TimePoint(), i.e. to have 0.0 time point)
 */
template<typename RealSpec>
class DiscreteTimeGrid {
public:
	/** Badge must be a single word! */
	static std::string badge() { return "DiscreteTimeGrid"; }

	typedef DiscreteTimeGrid Class;
	typedef DiscreteTimeGrid GridType;
	class TimePointType {
	public:
		typedef RealSpec RealType;

		/** Badge must be a single word! */
		static std::string badge() { return "DiscreteTimePoint"; }

		inline const GridType& grid() const { return this->m_grid; }
		explicit operator int() const { return this->toInt(); }
		inline int toInt() const { return this->m_i; }
		inline operator RealType() const { return *(m_grid.m_ptr_h) * RealType(m_i); }
		inline bool isZero() const { return m_i == 0 || m_grid.m_ptr_h == ptr_zero; }
		inline std::string show() { std::ostringstream oss; oss << (*m_ptr_h) << " * " << m_i; return oss.str(); }
		friend inline TimePointType operator+(TimePointType const & a, TimePointType const & b) {
			a.checkGridCompatible(b, "Cannot add time points");
			return TimePointType(a.isZero() ? b.m_grid : a.m_grid, a.m_i + b.m_i);
		}
		friend inline TimePointType operator-(TimePointType const & a, TimePointType const & b) {
			a.checkGridCompatible(b, "Cannot subtract time points");
			return TimePointType(a.isZero() ? b.m_grid : a.m_grid, a.m_i - b.m_i);
		}
		inline TimePointType operator+=(TimePointType const & b) {
			checkGridCompatible(b, "Cannot add time to this point");
			m_i += b.m_i;
			return *this;
		}
		inline TimePointType operator-=(TimePointType const & b) {
			checkGridCompatible(b, "Cannot subtract time from this point");
			m_i -= b.m_i;
			return *this;
		}
		/**
		 * Returns a point that is b steps after the point a on the grid.
		 * this is faster than adding two TimePoints, does not need to check compatibility!
		 */
		friend TimePointType operator+(TimePointType const & a, int b) {
			return TimePointType(a.m_grid, a.m_i + b);
		}
		/**
		 * Returns a point that is b steps before on the grid.
		 * this is faster than adding two TimePoints, does not need to check compatibility!
		 */
		friend inline TimePointType operator-(TimePointType const & a, int b) { return a + (-b); }
		/**
		 * Moves the point b steps further on the grid.
		 * this is faster than adding two TimePoints, does not need to check compatibility!
		 */
		inline TimePointType operator+=(int b) {
			m_i += b;
			return *this;
		}
		/**
		 * Moves the point is b steps back on the grid.
		 * this is faster than adding two TimePoints, does not need to check compatibility!
		 */
		inline TimePointType operator-=(int b) { return this->operator+=(-b); }
		/** this is getting -t from t */
		friend inline TimePointType operator-(const TimePointType& a) { return TimePointType(a.m_grid, - a.m_i); }
		friend inline TimePointType operator*(int const & a, TimePointType const & b) { return TimePointType(b.m_grid, b.m_i * a); }
		friend inline TimePointType operator*(TimePointType const & a, int const & b) { return b * a; }
		friend inline TimePointType operator*(RealType const & a, TimePointType const & b) { throw std::logic_error("DiscreteTimeGrid::TimePoints does not support multiplication by a real value! Only integer multiplication is allowed!"); }
		friend inline TimePointType operator*(TimePointType const & a, RealType const & b) { return b * a; }
		inline bool operator==(TimePointType const & other) const { return (&this->m_grid == &other.m_grid) && (this->m_i == other.m_i); }
		inline bool operator!=(TimePointType const & other) const { return !(*this == other); }
		TimePointType(TimePointType const & orig): m_grid(orig.m_grid), m_i(orig.m_i) {}
		TimePointType(): m_grid(GridType::trivialGrid), m_i(0) {}
		TimePointType& operator=(TimePointType const & orig) {
			this->checkGridCompatible(orig);
			if (m_grid.m_ptr_h == ptr_zero)
				throw std::logic_error("DiscreteTimeGrid::TimePointType: cannot assign to absolute zero-point");
			this->m_i = orig.m_i;
			return *this;
		}
		friend inline  bool operator<(TimePointType const & a, TimePointType const & b) {
			try { a.checkGridCompatible(b); } catch (std::logic_error& e) { return RealType(a) < RealType(b); }
			return a.m_i < b.m_i;
		}
		friend inline bool operator>(TimePointType const & a, TimePointType const & b) { a.checkGridCompatible(b); return b < a; }
		friend inline bool operator<=(TimePointType const & a, TimePointType const & b) { a.checkGridCompatible(b); return !(a > b); }
		friend inline bool operator>=(TimePointType const & a, TimePointType const & b) { a.checkGridCompatible(b); return !(b > a); }

		friend std::ostream& operator<<(std::ostream & out, TimePointType const & t) {
			out << RealType(t) << " := " << TimePointType::badge() << " ";
			out << t.m_i << " " << t.m_grid.h();
			return out;
		}
		friend std::istream& operator>>(std::istream & in, TimePointType & t) {
			RealType dump_repr; in >> dump_repr;
			helper_dump_badge(in); // dump ":="
			helper_dump_badge(in); // dump badge() = "DiscreteTimePoint"
			RealType tmp;
			in >> t.m_i >> tmp;
			if (tmp != *(t.m_grid.m_ptr_h)){ // check at least values
				std::ostringstream info;
				info << "TimePointType: Point from input is probably from some other grid. ";
				info << "Input h = " << tmp << " vs grid h = " << *(t.m_grid.m_ptr_h);
				throw std::logic_error(info.str());
			}
			return in;
		}
	protected:
		const GridType& m_grid;
		int m_i;
		void checkGridCompatible(TimePointType const & other, std::string extraInfo = "") const {
			if (this->isZero() || other.isZero()) return; 		// zero is compatible to all
			if (this->m_grid.m_ptr_h != other.m_grid.m_ptr_h){	// points must be physically the same
				std::ostringstream info;
				info << "DiscreteTimePoint: " << extraInfo << ". The points come from different grids (even if they have the same step size h)!";
				info << *(this->m_grid.m_ptr_h) << " vs " << *(other.m_grid.m_ptr_h) << "\n";
				info << "[" <<  this->m_grid.m_ptr_h << " vs " << other.m_grid.m_ptr_h << "]\n";
				throw std::logic_error(info.str());
			}
		}
		TimePointType(GridType const & grid, int i): m_grid(grid), m_i(i){}
		friend class DiscreteTimeGrid<RealType>;
	};
	friend class TimePointType;

	DiscreteTimeGrid(DiscreteTimeGrid const& orig): m_ptr_h(orig.m_ptr_h) {}
	DiscreteTimeGrid& operator=(DiscreteTimeGrid const& orig) { m_ptr_h = orig.m_ptr_h; return *this; }
	friend bool operator==(DiscreteTimeGrid const& a, DiscreteTimeGrid const& b) { return a.m_ptr_h == b.m_ptr_h; }
	friend bool operator!=(DiscreteTimeGrid const& a, DiscreteTimeGrid const& b) { return !(a == b); }
	DiscreteTimeGrid(): m_ptr_h(ptr_zero) {}
	DiscreteTimeGrid(const RealSpec& h): m_ptr_h(h == zero? ptr_zero : std::make_shared<RealSpec>(h)) {} // this way we always point out to the same 0!
	TimePointType point(int i) const { return TimePointType(*this, i); }
	TimePointType operator()(int i) const { return TimePointType(*this, i); }
	void split(RealSpec t, TimePointType& ti, RealSpec& epsi) const {
		int i = 0;
		if (m_ptr_h != ptr_zero){ i = closestSmallerInt(t / *m_ptr_h); }
		ti = TimePointType(this, i);
		epsi = t - RealSpec(ti);
	}
	RealSpec h() const { return *m_ptr_h; }
	static const RealSpec zero;  						// this way there is always exactly one physical zero we can use
	static const std::shared_ptr<RealSpec> ptr_zero;  	// used to initialize shared ptrs
	static const GridType trivialGrid;  				// this way there is always at least one grid
protected:
	std::shared_ptr<RealSpec> m_ptr_h;
	DiscreteTimeGrid(const std::shared_ptr<RealSpec> ptr_h): m_ptr_h(ptr_h) {} 	 // making it private, so no one tries to override
};
template<typename RealSpec>
const RealSpec DiscreteTimeGrid<RealSpec>::zero = 0;
template<typename RealSpec>
const std::shared_ptr<RealSpec>  DiscreteTimeGrid<RealSpec>::ptr_zero = std::make_shared<RealSpec>(DiscreteTimeGrid<RealSpec>::zero);
template<typename RealSpec>
const DiscreteTimeGrid<RealSpec> DiscreteTimeGrid<RealSpec>::trivialGrid(ptr_zero);

/**
 * computes Taylor sum given step and coefficients:
 * sum_{k=0}^{n} a[k] * step^k
 * a is an forward iterator (by copy, so we can modify it)
 * This iterator returs elements compatible with OutSpec (Scalar, Vector, etc)
 * OutSpec - scalar or vector to store result. 
 *
 * NOTE: out variable should be set to appropriate 0 outside of function!
 *
 * This is a naive implementation, that does not assume 
 * Taylor collection can be viewed from the back. See sumTaylorBackward()
 * for a more robust version with Horner's sumation algorithm.
 *
 * TODO: (NOT URGENT) there should be a version of this that accepts begin and end iterators.
 * TODO: (NOT URGENT) propose to implement this somewhere in CAPD itself
 * TODO: (NOT URGENT) (I have seen a lot of reimplementations there, DRY)
 */
template<typename ForwardIteratorInSpec, typename StepScalarSpec, typename OutSpec>
void sumTaylorForward(ForwardIteratorInSpec a, int n, StepScalarSpec step, OutSpec& out){
	StepScalarSpec powr = StepScalarSpec(1.0);	// this is current step^k
	for (int j = 0; j <= n; j++){				// sum all coefficents
		out += *a * powr;						// sum next element
		powr *= step;							//
		++a;
	}
} // sumTaylorForward()

/**
 * computes Taylor sum given step and coefficients:
 * sum_{k=0}^{n} a[k] * step^k
 * a is backward iterator, so we can take the advantage of Horner's algorithm
 * OutSpec - scalar or vector to store result. 
 *
 * NOTE: out variable should be set to appropriate 0 outside of function!
 *
 * This is a naive implementation, that does not assume 
 * Taylor collection can be viewed from the back. See sumTaylorBackward()
 * for a more robust version with Horner's sumation algorithm.
 *
 * TODO: (NOT URGENT) there should be a version of this that accepts begin and end iterators.
 * TODO: (NOT URGENT) propose to implement this somewhere in CAPD itself (I have seen a lot of reimplementations there, DRY)
 * TODO: (NOT URGENT) I myself does not use this to DRY... :(
 */
template<typename BackwardIteratorInSpec, typename StepScalarSpec, typename OutSpec>
void sumTaylorBackward(BackwardIteratorInSpec a, int n, StepScalarSpec step, OutSpec& out){			
	for (int j = 0; j <= n; j++){						
		out = *a + (step * out);
		--a;
	}
} // sumTaylorBackward()


template<typename MatrixType, typename size_type>
std::vector<MatrixType> extractDiagonalBlocks(MatrixType const& M, size_type d, size_type& offBlockDiagonalCount){
	std::vector<MatrixType> result;
	size_type I = 0;
	MatrixType B(d, d);
	size_type cols = M.numberOfColumns();
	size_type rows = M.numberOfRows();
	if (cols % d != 0 || rows % d != 0){
		std::ostringstream info;
		info << "extractDiagonalBlocks(): rows " << rows << " or cols " << cols << " are not multiple of block d " << d;
		throw std::range_error(info.str());
	}
	while (I < rows){
		for (size_type i = 0; i < d; i++){
			for (size_type J = 0; J < cols; J++){
				if (I <= J && J < I + d)
					B[i][J-I] = M[I+i][J];
				else if (M[I+i][J] != 0.)
					++offBlockDiagonalCount;
			}
		}
		I += d;
		result.push_back(B);
	}
	return result;
}


} // namespace ddes
} // namespace capd


#endif /* _CAPD_DDES_DDECOMMON_H_ */
