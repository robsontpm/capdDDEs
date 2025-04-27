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

#ifndef _CAPD_DDES_NONRIGOROUSTAYLORSOLVER_H_
#define _CAPD_DDES_NONRIGOROUSTAYLORSOLVER_H_

#include <stdexcept>
#include <string>
#include <capd/ddes/DDECommon.h>
#include <capd/ddes/FunctionalMap.h>

namespace capd{
namespace ddes{

/**
 * a class to hold additional data to compute Jacobian of the
 * (nonrigorous / approximate) flow w.r.t. initial data.
 * Used in the Taylor integrator to do the C^1 computations.
 *
 * It has the same interface as the VectorSpec, but holds extra data.
 * The current implementation holds the matrix D of shape d x M
 * as a collection of blocks d x d. This is beneficial from the
 * implementation point of view in the Taylor integration
 * and in the extraction of data. Basically, each coefficient coeff
 * in Taylor representation is a d dimensional vector, so the
 * derivative of coeff w.r.t. other coeff is a d x d matrix. The derivative
 * of a coeff w.r.t. collection of coeffs is a collection of such matrices.
 */
template<typename VectorSpec, typename MatrixSpec>
class VectorWithJacData : public VectorSpec {
public:
	typedef MatrixSpec MatrixType;
	typedef typename MatrixSpec::ScalarType ScalarType;
	typedef VectorWithJacData<VectorSpec, MatrixSpec> Class;
	typedef VectorSpec BaseClass;
	typedef typename BaseClass::size_type size_type;
	typedef std::vector<MatrixSpec> MatrixStorageType;

	// DEV NOTE: I deliberately left the implementations here, they are very short
	// DEV NOTE: and i like to have them in one place. Please do not refactor this.

	VectorWithJacData(): BaseClass() {}
	VectorWithJacData(size_type d): BaseClass(d) {}
	VectorWithJacData(const VectorWithJacData& orig): BaseClass(orig), m_Jac(orig.m_Jac) {}
	VectorWithJacData(const VectorSpec& orig): BaseClass(orig) {}
	VectorWithJacData(const VectorSpec& v, const std::vector<MatrixSpec>& D): BaseClass(v), m_Jac(D) { }
	VectorWithJacData(const VectorSpec& v, const MatrixSpec& D): BaseClass(v) {
		try {
			setMatrix(D);
		} catch (std::logic_error& e){
			throw rethrow("VectorWithJacData::__construct__(): ", e);
		}
	}

	/**
	 * sets the data with a given matrix,
	 * matrix must be d x M, where d is the dimension of
	 * this vector. Otherwise std::logic_error is thrown.
	 */
	VectorWithJacData& setMatrix(MatrixType const& D){
		size_type d = this->dimension();
		if (D.numberOfRows() != d){
			std::ostringstream info;
			info << "VectorWithJacData::setMatrix(): D has bad number of rows, ";
			info << "expected: " << d << ", is: " << D.numberOfRows();
			throw std::logic_error(info.str());
		}
		size_type I = 0;
		m_Jac.clear();
		while(I < D.numberOfColumns()){
			MatrixSpec M(d, d);
			// first go along columns (i, I),
			// see that I changes same number as i, so after loop I is increased by d
			// then we go with j along each row for each column.
			// this way we partitions matrix D into matrices M of dim (d,d)
			for (size_type i = 0; i < d && I < D.numberOfColumns(); ++i, ++I)
				for (size_type j = 0; j < d; j++)
					M[j][i] = D[j][I];
			m_Jac.push_back(M);
		}
		return *this;
	}

	/** standard assign */
	VectorWithJacData& operator=(const VectorWithJacData& orig){
		BaseClass::operator=(BaseClass(orig));
		m_Jac = orig.m_Jac;
		return *this;
	}

	/** standard assign from the base class, with empty matirx data. */
	VectorWithJacData& operator=(const BaseClass& orig){
		BaseClass::operator=(orig);
		m_Jac.clear();
		return *this;
	}

	/**
	 * standard multiplication by a scalar. It also modifies Matrix data D,
	 * as if D is the derivative of x, then cD is derivative of cx.
	 */
	VectorWithJacData& operator*=(ScalarType const& c){
		this->VectorSpec::operator*=(c);
		for (auto i = m_Jac.begin(); i != m_Jac.end(); ++i) (*i) *= c;
		return *this;
	}

	/**
	 * adds two vectors with data.
	 * WARNING: the matrix data must be compatible!
	 *          If not, there will be probably an exception thrown
	 *          by the CAPD internal implementation.
	 */
	VectorWithJacData& operator+=(VectorWithJacData const& other){
		this->VectorSpec::operator+=(other);
		auto i = m_Jac.begin();
		auto o = other.m_Jac.begin();
		for (; i != m_Jac.end() && o != other.m_Jac.end(); ++i, ++o) (*i) += (*o);
		return *this;
	}

	/** conversion to the basic vector type */
	operator VectorSpec() { return *this; }
	/** converts the sequence of Matrices d x d, into one big matrix d x M */
	operator MatrixSpec() {
		size_type d = this->dimension();
		size_type N0 = d * m_Jac.size();
		size_type I = 0;
		MatrixSpec M(d, N0);
		for (auto jaci = m_Jac.begin(); jaci != m_Jac.end(); ++jaci)
			for (size_type i = 0; i < d; ++i, ++I)
				for (size_type j = 0; j < d; ++j)
					M[j][I] = (*jaci)[j][i];
		return M;
	}

	/** returns the internal matrix data for a protected variable */
	const MatrixStorageType& getMatrixData() const { return m_Jac; }
	/** returns the internal matrix data and allows modyfication. */
	MatrixStorageType& getMatrixData() { return m_Jac; }

protected:
	std::vector<MatrixSpec> m_Jac;
};

/**
 * This is nonrigorous Taylor method for general FDE/DDEs,
 * such that they can provide some basic Jet data on their
 * representation (@see computeDDECoefficients in an exemplary class
 * from @see BasicDiscreteDelaysFunctionalMap.h).
 *
 * This is the implementation of an algorithm given
 * in a series of FoCM papers (2018 and 2024). See them for
 * mathematical background.
 */
template<typename FunctionalMapSpec>
class DDENonrigorousTaylorSolver {
public:
	/** naming required by CAPD and other ddes codes. It simplifies coding */
	typedef FunctionalMapSpec FunctionalMapType;
	/**
	 * naming required by CAPD and other ddes codes. It simplifies coding
	 * this is not the Vector Field per se, but it should suffice for now, and I think it is what is called in DDEs literature.
	 * for now, it is just to have some compatibility with CAPD if we plan to include this library in CAPD
	 */
	typedef FunctionalMapType VectorFieldType;
	/** ddes compatible. Naming required by CAPD and other ddes codes. It simplifies coding */
	typedef typename FunctionalMapType::CurveType SolutionCurveType; // (For now, I will stick to my version, I think CAPD has here a problem with naming convention)
	/** CAPD compatible. Naming required by CAPD and other ddes codes. It simplifies coding */
	typedef typename FunctionalMapType::CurveType SolutionCurve;
	/** CAPD compatible. Naming required by CAPD and other ddes codes. It simplifies coding */
	typedef typename FunctionalMapType::CurveType CurveType;
	/** CAPD compatible. Naming required by CAPD and other ddes codes. It simplifies coding */
	typedef typename FunctionalMapType::JetType JetType;
	/** CAPD compatible. Naming required by CAPD and other ddes codes. It simplifies coding */
	typedef typename FunctionalMapType::DataType DataType; 	// it is in fact SetType
	/** CAPD compatible. Naming required by CAPD and other ddes codes. It simplifies coding */
	typedef typename FunctionalMapType::VectorType VectorType;
	/** CAPD compatible. Naming required by CAPD and other ddes codes. It simplifies coding */
	typedef typename FunctionalMapType::MatrixType MatrixType;
	/** Naming required by CAPD and other ddes codes. It simplifies coding */
	typedef typename FunctionalMapType::TimePointType TimePointType;
	/** Naming required by CAPD and other ddes codes. It simplifies coding */
	typedef typename FunctionalMapType::RealType RealType;
	/** CAPD compatible. Naming required by CAPD and other ddes codes. It simplifies coding */
	typedef typename FunctionalMapType::ScalarType ScalarType;
	/** Naming required by CAPD and other ddes codes. It simplifies coding */
	typedef typename FunctionalMapType::size_type size_type;
	/** This is needed by C^1 computations only */
	typedef typename FunctionalMapType::VariableStorageType VariableStorageType;
	/** This is needed by C^1 computations only */
	typedef typename FunctionalMapType::JacobianStorageType JacobianStorageType;
	/** This is needed by C^1 computations only */
	typedef typename FunctionalMapType::ValueStorageType ValueStorageType;

	/**
	 * basic constructor.
	 * We assume that no solver can be created without giving the FDE/DDE
	 * in form of FunctionalMapType (@see BasicDiscreteDelaysFunctionalMap.h
	 * for an exemplary class that can be used here).
	 */
	DDENonrigorousTaylorSolver(FunctionalMapType const& map, size_type maxOrder = 20): m_map(map), m_maxOrder(maxOrder) {}
	/** copy constructor. NOTE: there is also a default assign operator=(). */
	DDENonrigorousTaylorSolver(DDENonrigorousTaylorSolver const& solver): m_map(solver.m_map), m_maxOrder(solver.m_maxOrder) {}


	/**
	 * Basic form of the operator to do one step. It has additional parameter
	 * that allows to do many steps at once (basically computing a time map, but only for the full steps
	 * of size h from the grid).
	 * By default it only do one step, so it computes approximate solution on [t_0, t_0 + h]
	 *
	 * @param in_out_curve
	 * 		this is the representation of a solution on some interval $[-t_past, t_0]$,
	 * 		after the procedure, it will contain the representation over $[-t_past, t_0 + steps*h]$
	 * @param steps number of full steps of size h to do, default is just 1
	 *
	 * NOTE: You can supply TimePointType as steps (see other function)!
	 * So the following code is valid:
	 *      int p = 128;			 // 128 points per one delay
	 *      Grid grid(1.0 / p);		 // make a grid with step 2 / p (we assume tau = 1 in the equation)
	 *		auto tau = grid(p); 	 // get a point from the grid that corresponds to the delay
	 *		solver(solution, tau*4); // integrate on [t_0, t_0 + 4\tau]
	 */
	void operator()(CurveType& in_out_curve, size_type steps=1);

	/*
	 * This allows to supply TimePointType as the destination time!
	 *
	 * @param in_out_curve
	 * 		this is the representation of a solution on some interval $[-t_past, t_0]$,
	 * 		after the procedure, it will contain the representation over $[-t_past, t_0 + T]$
	 * @param T the time to extend the solution
	 * So the following code is valid:
	 *      int p = 128;			 // 128 points per one delay
	 *      Grid grid(1.0 / p);		 // make a grid with step 2 / p (we assume tau = 1 in the equation)
	 *		auto tau = grid(p); 	 // get a point from the grid that corresponds to the delay
	 *		solver(solution, tau*4); // integrate on [t_0, t_0 + 4\tau]
	 */
	inline void operator()(CurveType& in_out_curve, TimePointType const&T){
		this->operator()(in_out_curve, int(T));
	}

	/**
	 * This does one step of integration with computing the variational
	 * equation. It is used to get the so called C^1 computations,
	 * i.e. computing the derivative of the flow w.r.t. initial data.
	 * Note this is rather cumbersome to use by the End User, so
	 * we propose to use the specialized class to get the derivative of the
	 * solution w.r.t. initial data, such as DDEBasicPoincareMap or
	 * DDETimeMap.
	 *
	 * NOTE: (important!)
	 * 		we assume that DataType of Jet can hold value and Matrix.
	 *      We provide candidate for this: VectorWithJacobianData
	 *      you cannot use this function with the basic data structure for SolutionCurve
	 *      this is also motivated with the optimization in computation of rhs
	 *      (if we do not use all past data, but only e.g. single discrete delay)
	 *      the VectorWithJacobianData is to bind value of a variable with the
	 *      derivative of coefficients w.r.t. initial data.
	 *      Please check examples on how to use this in your code.
	 */
	void operator()(CurveType& in_out_curve, JacobianStorageType& D);

	/** returns the functional map defining the r.h.s. of the equation (with delays and everything) */
	FunctionalMapType const& getMap() const { return m_map; }
	/** returns the functional map defining the r.h.s. of the equation (with delays and everything) */
	FunctionalMapType & getMap() { return m_map; }

	/**
	 * This extracts the variational matrix
	 * (i.e. the derivative of the solution with respect to initial segment)
	 * computed in the curve during the integration procedure.
	 *
	 * NOTE: the solution can be unaware about having "C^1-data" needed to obtain
	 * the variational matrix. Therefore, the procedure is not part of the
	 * solution class, but it is tied to the Solver class, which produces
	 * along the way the C^1 data if needed.
	 *
	 * See the FoCM papers for more details. But basically,
	 * here we treat the curve as a solution depending on the initial
	 * data curve = curve(initial). We want compute \partial{curve}\over\partial{u}(initial),
	 * where $u$ is from a functional space, so we have a Frechet derivative problem
	 * and hard to treat explicitly. But on the other hand, we know that
	 * u \equiv x \in \R^M, that is, there is a vector x in finite space
	 * that represents the u. The same for curve(u), it is represented by y.
	 * Therefore we can work on y(x), a map from \R^M to \R^{M'} (it is possible that)
	 * $M' \ne M$. And me can compute a matrix of the shape M' x M that
	 * decodes the partial derivatives \partial y_i \over \partial x_j.
	 * And this is what fullV below is.
	 *
	 * The reducedV is a technical thing, in fact, I use it only when computing
	 * M x M matrix, that decodes subset of the fullV that coresponds to the
	 * projection of y(x) to the space of representation of x.
	 *
	 * NOTE: the end user should hardly need to call this procedure
	 *
	 * For more information on the possible use:
	 * @see extractVariationalMatrix(SolutionCurve const&, MatrixType&, MatrixType&) const
	 */
	void extractVariationalMatrix(
	            SolutionCurve const& curve,
	            MatrixType& fullV, MatrixType& reducedV,
	            std::vector<size_type> reducedShape, int p_howFar = -1) const;

	/**
	 * This extracts the variational matrix from a given curve. This is
	 * the basic version, that extracts all data for all coefficients
	 * over the whole curve. It might be good for images of Poincare or TIme maps,
	 * but for a whole solution it might give you some weirdly shaped matrix
	 * without any meaning. You can get more control in the other version.
	 *
	 * For example (this is pseudocode, you need to setup computations first):
	 * Grid grid(1./100.);
	 * DDETimeMap T(...);
	 * C1SolutionCurve X(...);
	 * capd::DMatrix V;
	 * auto PX = T(grid(100), X); // move by time t = 1 = 100*grid.h()
	 * P.getSolver().extractVariationalMatrix(X, V);
	 *
	 * In the case of PoincareMaps it is better to use:
	 * capd::DMatrix DP;
	 * DDEBasicPoincareMap P(...);
	 * auto PX = P(X, DP);
	 * which return the true Jacobian of the map, corrected by the terms related to the return time t_p
	 * See the FoCM papers for more details.
	 *
	 * For more information:
	 * @see extractVariationalMatrix(SolutionCurve const&, MatrixType&, MatrixType&, std::vector<size_type>, int) const
	 */
	void extractVariationalMatrix(SolutionCurve const& curve, MatrixType& fullV) const;

protected:
	/**
	 * same as the other oneStep, but without computing JacPhi
	 * NOTE: out_Phi_z is a collection, for possible future use,
	 *       I assume that there might be more than only value of x(t+h),
	 *       maybe from different algorithms and then we will combine them?
	 *       Or maybe we will compute more coefficients? We will see.
	 *
	 * It is designed to be used internally. Operators() and epsilonShift() use this.
	 *
	 * @param in_t0 the current time
	 * @param in_th next time, should be that th \approx t0 + h, note we want to have this as this is from grid, and t0 + h might be not (approximation)!
	 * @param in_h  step size, see above
	 * @param in_u  this is a representation of the past history that is needed to compute one step. Usually m_dynsys takes care of providing this parameter by collectComputationData()
	 * @param out_Phi_coeffs_t0 here will be stored the coefficients of the solution at t = t0
	 * @param out_Phi_z here will be stored the value(s) of the representation at t = th, usually just one = the value of solution = x(th)
	 */
	void oneStep(
			TimePointType const& 		in_t0,
			TimePointType const&		in_th,
			RealType const&				in_h,
			ValueStorageType const&		in_u,
			ValueStorageType& 			out_Phi_coeffs_t0,
			ValueStorageType& 			out_Phi_z);

	/**
	 * this function produces approximation to the solution
	 * after time h = distance between grid points in the input curve.
	 * The var out_Phi_coeffs_t0 will hold the Jet coefficients at point
	 * t0, then out_JacPhi_coeffs_t0 will be the Jacobian of
	 * coeffs w.r.t. variables in out_u.
	 * Then out_Phi_x[0] is the value at t = t0+h computed by the Taylor method,
	 * i.e. evaluation of out_Phi_coeffs_t0(t0+h). out_JacPhi_x[0] is
	 * like in Taylor method for ODEs.
	 *
	 * It is designed to be used internally. Operators() and epsilonShift() use this.
	 *
	 * NOTE: out_Phi_z is a collection, see other version of oneStep for more info.
	 * NOTE: in out_JacPhi_z I need storage for as many as there will be entries in out_Phi_z
	 *
	 * @param out_JacPhi_coeffs_t0
	 * 		this is a collection of matrices o shape d x M
	 *		(M = size of the representation of the initial function, see papers)
	 *		it stores the derivative \partial{coef}\over\partial{initial}(initial),
	 *		i.e. the derivative of the given coefficient w.r.t. initial data.
	 *		It looks terrible written like this, but see FoCM 2024 paper for proper
	 *		mathematical treatment.
	 * @param out_JacPhi_z
	 * 		similar as the above, but the derivative of the value at th w.r.t. initial data.
	 *
	 * NOTE: See other version of this procedure to get info on other parameters.
	 */
	void oneStep(
			TimePointType const& 		in_t0,
			TimePointType const&		in_th,
			RealType const&				in_h,
			ValueStorageType const&		in_v,
			VariableStorageType const&	in_u,
			JacobianStorageType const&	in_dvdu,
			ValueStorageType& 			out_Phi_coeffs_t0,
			JacobianStorageType& 		out_JacPhi_coeffs_t0,
			ValueStorageType& 			out_Phi_z,
			JacobianStorageType& 		out_JacPhi_z);

private:
	/**
	 * definition of the actual FDE/DDE, with delays ant everything.
	 * @see BasicDiscreteDelaysFunctionalMap.h for the representative
	 * of such a class.
	 */
	FunctionalMapType m_map;
	/**
	 * maximal order of jet produced at each step.
	 * This might be reduced by the max allowed order due to the Curve continuity.
	 * See FoCM papers for more information.
	 */
	size_type m_maxOrder;
};



} // namespace ddes
} // namespace capd


#endif /* _CAPD_DDES_NONRIGOROUSTAYLORSOLVER_H_ */
