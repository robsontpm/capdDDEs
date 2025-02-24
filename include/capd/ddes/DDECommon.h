/*
 * This file constitutes part of DDEs rigorous integration framework developed
 * in PhD Thesis under supervision of prof. Piotr Zgliczynski:
 *
 * 		"Rigorous Integration of Delay Differential Equations", Jagiellonian University, 2015
 *
 * When using in scientific work please consider citing my articles (preferred),
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

/** this is used in many places. Pretty technical. Nothing to do with DDEs actually. */
template<typename PtrT>
void helper_safe_delete(PtrT& item, bool is_owner){
	if (item && is_owner){
		delete item;
		item = 0;
	}
}

/** this is used in many places. Pretty technical. Nothing to do with DDEs actually. */
template<typename PtrT>
void helper_safe_array_delete(PtrT& item, bool is_owner){
	if (item && is_owner){
		delete[] item;
		item = 0;
	}
}

/** this is used in many places. Pretty technical. Nothing to do with DDEs actually. */
void helper_dump_line(std::istream & in);

/** this is used in many places. Pretty technical. Nothing to do with DDEs actually. */
void helper_dump_badge(std::istream & in);

/** helper function, to be used to simply override versions for double (numerical nonrigorous) and Interval (rigorous) */
template<typename RealSpec>
RealSpec ecloseStep(RealSpec h){ return RealSpec(0.0, 1.0) * h; }
/** nonrigorous version for double numerical only method */
double ecloseStep(double h);

/** helper function for pretty printing. */
template<typename RealSpec>
std::string showEnclosedInterval(RealSpec a, RealSpec b){ std::ostringstream oss; oss << "[" << a.leftBound() << ", " << b.rightBound() << ")"; return oss.str(); }
/** nonrigorous version for double numerical only method */
std::string showEnclosedInterval(double a, double b);

/**
 * this is used in many places. Pretty technical. Nothing to do with DDEs actually.
 * It allows to have some kind of stack trace as in Python, but in C++. Now
 * probably there is some library for this in boost/std.
 *
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
 * this is used in many places. Pretty technical. Nothing to do with DDEs actually.
 * It allows to have some kind of stack trace as in Python, but in C++. Now
 * probably there is some library for this in boost/std.
 *
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

/**
 * This represents time points on a grid i * h, i - integer, h - real
 * The mechanism is that, the grid holds a pointer to a scalar (physical quantity)
 * which might be an interval, but our interpretation is that: the scalar
 * represents exact point on the grid. If another point comes from the
 * grid then it is also an exact point on the grid, even if we have only estimation on it.
 * We prevent two grids with different physical step sizes h to mix, e.g. adding together, or
 * going into the same piecewise Taylor curve representation. By evaluating the curve
 * on the grid points, we can warrant some other properties or use better estimates.
 *
 * In general, you create one grid for all your computations and use this grid to
 * supply all time points for all other objects. You can use
 * capd::ddeshelper::(Non)RigorousHelper to handle this.
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
		/** returns the reference to a grid this points belongs to */
		inline const GridType& grid() const { return this->m_grid; }
		/** returns the number of the point in the grid, i.e. i from t_i = ih */
		explicit operator int() const { return this->toInt(); }
		/** returns the number of the point in the grid, i.e. i from t_i = ih */
		inline int toInt() const { return this->m_i; }
		/** conversion operator. It gives the evaluation of the point as a basic scalar type, i.e. the value of t_i = ih */
		inline operator RealType() const { return *(m_grid.m_ptr_h) * RealType(m_i); }
		/** tests if this point is a true zero */
		inline bool isZero() const { return m_i == 0 || m_grid.m_ptr_h == ptr_zero; }
		/** shows the human readable representation, without badge */
		inline std::string show() { std::ostringstream oss; oss << *(m_grid.m_ptr_h) << " * " << m_i; return oss.str(); }
		/** add two time points together. It tests if they come from the same grid */
		friend inline TimePointType operator+(TimePointType const & a, TimePointType const & b) {
			a.checkGridCompatible(b, "Cannot add time points");
			return TimePointType(a.isZero() ? b.m_grid : a.m_grid, a.m_i + b.m_i);
		}
		/** difference between two time points. It tests if they come from the same grid */
		friend inline TimePointType operator-(TimePointType const & a, TimePointType const & b) {
			a.checkGridCompatible(b, "Cannot subtract time points");
			return TimePointType(a.isZero() ? b.m_grid : a.m_grid, a.m_i - b.m_i);
		}
		/** add the value of a time point to this one. It tests if they come from the same grid */
		inline TimePointType operator+=(TimePointType const & b) {
			checkGridCompatible(b, "Cannot add time to this point");
			m_i += b.m_i;
			return *this;
		}
		/** substracts the value of a time point to this one. It tests if they come from the same grid */
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
		/** Moves the point to the next point in the grid (preincrement), faster than increment! */
		inline const TimePointType& operator++() { ++m_i; return *this; }
		/** Moves the point to the next point in the grid (postincrement), slower than predecrement! */
		inline TimePointType operator++(int) { TimePointType cpy(*this); ++m_i; return cpy; }
		/** Moves the point to the previous point in the grid (predecrement), faster than increment! */
		inline const TimePointType& operator--() { --m_i; return *this; }
		/** Moves the point to the next point in the grid (postdecrement), slower than predecrement! */
		inline TimePointType operator--(int) { TimePointType cpy(*this); --m_i; return cpy; }
		/**
		 * Moves the point is b steps back on the grid.
		 * this is faster than adding two TimePoints, does not need to check compatibility!
		 */
		inline TimePointType operator-=(int b) { return this->operator+=(-b); }
		/** this is getting -t from t */
		friend inline TimePointType operator-(const TimePointType& a) { return TimePointType(a.m_grid, - a.m_i); }
		/** we only support multiplication by integers */
		friend inline TimePointType operator*(int const & a, TimePointType const & b) { return TimePointType(b.m_grid, b.m_i * a); }
		/** we only support multiplication by integers */
		friend inline TimePointType operator*(TimePointType const & a, int const & b) { return b * a; }
		/** we only support multiplication by integers */
		friend inline TimePointType operator*(RealType const & a, TimePointType const & b) { throw std::logic_error("DiscreteTimeGrid::TimePoints does not support multiplication by a real value! Only integer multiplication is allowed!"); }
		/** we only support multiplication by integers */
		friend inline TimePointType operator*(TimePointType const & a, RealType const & b) { return b * a; }
		/** NOTE: If the grid step is not the same physical value, than they are different. */
		inline bool operator==(TimePointType const & other) const { return (this->m_grid.m_ptr_h == other.m_grid.m_ptr_h) && (this->m_i == other.m_i); }
		/** NOTE: If the grid step is not the same physical value, than they are different. */
		inline bool operator!=(TimePointType const & other) const { return !(*this == other); }
		/** copy constructor */
		TimePointType(TimePointType const & orig): m_grid(orig.m_grid), m_i(orig.m_i) {}
		/** this creates a dummy time point on a grid that contains only point 0. You should avoid this, as you cannot change this point to anything else! */
		TimePointType(): m_grid(GridType::trivialGrid), m_i(0) {}
		/** this copies the point */
		TimePointType& operator=(TimePointType const & orig) {
			this->checkGridCompatible(orig);
			if (m_grid.m_ptr_h == ptr_zero)
				throw std::logic_error("DiscreteTimeGrid::TimePointType: cannot assign to absolute zero-point");
			this->m_i = orig.m_i;
			return *this;
		}
		/** works only if the grids are compatible, i.e. have the same physical step h. Otherwise throws std::logic_error. */
		friend inline  bool operator<(TimePointType const & a, TimePointType const & b) {
			try { a.checkGridCompatible(b); } catch (std::logic_error& e) { return RealType(a) < RealType(b); }
			return a.m_i < b.m_i;
		}
		/** works only if the grids are compatible, i.e. have the same physical step h. Otherwise throws std::logic_error. */
		friend inline bool operator>(TimePointType const & a, TimePointType const & b) { a.checkGridCompatible(b); return b < a; }
		/** works only if the grids are compatible, i.e. have the same physical step h. Otherwise throws std::logic_error. */
		friend inline bool operator<=(TimePointType const & a, TimePointType const & b) { a.checkGridCompatible(b); return !(a > b); }
		/** works only if the grids are compatible, i.e. have the same physical step h. Otherwise throws std::logic_error. */
		friend inline bool operator>=(TimePointType const & a, TimePointType const & b) { a.checkGridCompatible(b); return !(b > a); }

		/** standard output, it is compatible with standard input */
		friend std::ostream& operator<<(std::ostream & out, TimePointType const & t) {
			out << RealType(t) << " := " << TimePointType::badge() << " ";
			out << t.m_i << " " << t.m_grid.h();
			return out;
		}
		/** standard input, it is compatible with standard output */
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
		inline bool sameGrid(TimePointType const& other){ return this->grid() == other.grid(); }
	protected:
		const GridType& m_grid;
		int m_i;
		/** points are compatible only if one of them is 0, or if the grids have the same physical h */
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
		/** this is protected to prevent user of creating TimePoints without defining grid first. */
		TimePointType(GridType const & grid, int i): m_grid(grid), m_i(i){}
		friend class DiscreteTimeGrid<RealType>;
	};
	friend class TimePointType;

	/** copy constructor */
	DiscreteTimeGrid(DiscreteTimeGrid const& orig): m_ptr_h(orig.m_ptr_h) {}
	/** copy operator */
	DiscreteTimeGrid& operator=(DiscreteTimeGrid const& orig) { m_ptr_h = orig.m_ptr_h; return *this; }
	/** two grids are equal if they have the same physical h */
	friend bool operator==(DiscreteTimeGrid const& a, DiscreteTimeGrid const& b) { return a.m_ptr_h == b.m_ptr_h; }
	/** two grids are equal if they have the same physical h */
	friend bool operator!=(DiscreteTimeGrid const& a, DiscreteTimeGrid const& b) { return !(a == b); }
	/** this creates a dummy grid that contain the only point 0, with step size 0 */
	DiscreteTimeGrid(): m_ptr_h(ptr_zero) {}
	/** the most common contructor, it constructs a grid with step size h */
	DiscreteTimeGrid(const RealSpec& h): m_ptr_h(h == zero? ptr_zero : std::make_shared<RealSpec>(h)) {} // this way we always point out to the same 0!
	/** returns a new time point on this grid t_i = h*i */
	TimePointType point(int i) const { return TimePointType(*this, i); }
	/** returns a new time point on this grid t_i = h*i */
	TimePointType operator()(int i) const { return TimePointType(*this, i); }
	/** checks if a given point is from that grid */
	bool has(TimePointType const& ti) const { return *this == ti.grid(); }
	/**
	 * it tries to compute the representation of t as $t_i + \epsilon$,
	 * where $t_i = i * h$. It works only if $t_i$ is from this grid.
	 */
	void split(RealSpec t, TimePointType& ti, RealSpec& epsi) const {
		int i = 0;
		if (m_ptr_h != ptr_zero){ i = closestSmallerInt(t / *m_ptr_h); }
		ti = TimePointType(*this, i);
		if (RealSpec(ti + 1) <= t) ++ti;
		epsi = t - RealSpec(ti);
	}
	/** returns the value of step size h */
	RealSpec h() const { return *m_ptr_h; }
	static const RealSpec zero;  						///< this way there is always exactly one physical zero we can use
	static const std::shared_ptr<RealSpec> ptr_zero;  	///< used to initialize shared ptrs
	static const GridType trivialGrid;  				///< this way there is always at least one grid
protected:
	/** this is the pointer to a physical copy. The pointer is used to be able to have a physical reference and also to be able to modify it if needed */
	std::shared_ptr<RealSpec> m_ptr_h;
	/** making it private, so no one tries to override my check if the h is zero form the other constructor. */
	DiscreteTimeGrid(const std::shared_ptr<RealSpec> ptr_h): m_ptr_h(ptr_h) {}
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
 * a is backward iterator (that supports operator--()!),
 * so we can take the advantage of Horner's algorithm
 * OutSpec - scalar or vector to store result. 
 *
 * NOTE: out variable should be set to appropriate 0 outside of function!
 *
 * This is a naive implementation, that does not assume 
 * Taylor collection can be viewed from the back. See sumTaylorBackward()
 * for a more robust version with Horner's sumation algorithm.
 *
 * TODO: (NOT URGENT) Name might be misleading, as somebody might think of using vector.rbegin(), but in that case the iterator goes back by ++, not -- and this procedure returns rubbish! Rethink. Or maybe just add some docs notes?
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
