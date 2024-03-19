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

#ifndef _CAPD_DDES_FUNCTIONALMAP_H_
#define _CAPD_DDES_FUNCTIONALMAP_H_

#include <capd/ddes/DDECommon.h>
#include <capd/basicalg/TypeTraits.h>
#include <stdexcept>
#include <string>

namespace capd{
namespace ddes{

/**
 * * TODO: rewrite DOCS
 * this is the interface of a Functional map, that can be evaluated
 * only by "querying" the curve about its jets at various points in the past.
 * it can produce jet at a given time t together with partial derivatives
 * w.r.t. to internal variables from curve used in the computation by
 * recurrent formula for a functional equation:
 *
 * 		x'(t) = f(t, x)
 *
 * where x is the solution curve in \R^d defined for times smaller than
 * t (up to some maximal past time usually).
 *
 * A good example of such a functional map is discrete delay case,
 * e.g. F(t, x) = f(t, x(t), x(t-tau)), f : \R \times \R^{2d}->\R
 * (easily extended to more delays). (Probably) the functionals
 * that include integral over the past: F(t, x) = f(t, \int_{t-tau}^{t} x(s) ds)
 * will also fall into this category.
 *
 * The concrete implementations will be available later.
 *
 * Implementation note: when implementing concrete class of
 * this interface, do not forget to put
 * 	 using BaseClass::operator();
 *	 using BaseClass::computeDDECoefficients;
 * in the public part of your class to gain access to default
 * implementations of some functions.
 *
 * TODO: (NOT URGENT, FUTURE, RETHINK) it would be better not to pass containers to computeDDECoefficients, but maybe iterators to given types? It could be used to compute in a given place, i.e. at an already created Jet?
 */
template<typename SolutionCurveSpec, typename JetSpec=typename SolutionCurveSpec::JetType>
class DDEBasicFunctionalMap {
public:
	typedef SolutionCurveSpec CurveType;
	typedef JetSpec JetType;
	typedef typename CurveType::RealType RealType;
	typedef typename CurveType::TimePointType TimePointType;
	typedef typename CurveType::MatrixType MatrixType;
	typedef typename CurveType::VectorType VectorType;
	typedef typename CurveType::ScalarType ScalarType;
	typedef typename CurveType::size_type size_type;
	typedef typename CurveType::DataType DataType;

	/**
	 * This type will store all variables from the input curve
	 * used to evaluate the map and to produce recurrently
	 * the jet at a given time. It will be stored in native
	 * SetType to allow usage of Lohner-type algorithms
	 * for controlling wrapping effect later.
	 * Each element of this collection would be a set
	 * in \R^d dimensional space.
	 */
	typedef std::vector< DataType > VariableStorageType;
	/**
	 * similar to VariableStorageType but does not care about
	 * the internal representation. It will be used to
	 * store the Jet coefficients only in output
	 * (we do not need to care about the time point and we
	 * can freely change order (i.e. size of the ValueStorageType)
	 */
	typedef std::vector< VectorType > ValueStorageType;
	/**
	 * This type will store partial derivatives with respect
	 * to variables from input curve used in the computations.
	 * Each matrix will be of M(d, d) type.
	 */
	typedef std::vector< std::vector<MatrixType> > JacobianStorageType;

	/** helper to convert from collection of sets to collection of vectors */
	static void convert(VariableStorageType const& u, ValueStorageType &v){
		v.clear();
		for (auto it = u.begin(); it != u.end(); ++it)
			v.push_back(VectorType(*it));
	}

	/**
	 * helper to convert from collection of sets to collection of vectors
	 * NOTE: name is changed, as ValueStorageType might be VariableStorageType
	 *       in some instances, and this would make compilation ambiguous.
	 */
	static void deconvert(ValueStorageType const& u, VariableStorageType &v){
		v.clear();
		for (auto it = u.begin(); it != u.end(); ++it)
			v.push_back(DataType(*it));
	}

	/** output dimension of the map. THIS FUNCTION NEEDS TO BE IMPLEMENTED IN DERIVED CLASSES. */
	virtual size_type imageDimension() const = 0;
	/**
	 * computes value of the map at a given time point for a given solution curve.
	 * THIS FUNCTION NEEDS TO BE IMPLEMENTED IN DERIVED CLASSES.
	 */
	virtual VectorType operator()(const TimePointType& t0, const CurveType& x) const = 0;

	/**
	 * Collects all data from curve necessary for the computation of the value.
	 * This does not collect ENCLOSURE DATA required in rigorous computations.
	 *
	 * We use  collection of DataType items as the output, we ASSUME that DataType elements are
	 * used to describe function by CurveType. Therefore, in the output, we could
	 * relate how inputs (VariableStorage) corresponds to outputs (ValueStorage)
	 * mainly by Jacobian of the map at Variables (JacobianStorageType).
	 *
	 * Then, the main function that needs to be implemented in the Map
	 * are those depending on the VariableStorageType. Other can have
	 * default implementations and are for users convenience.
	 *
	 * THIS FUNCTION NEEDS TO BE IMPLEMENTED IN DERIVED CLASSES.
	 *
	 * Note: dt might seem redundant here, but for some reasons it might be helpful, i.e. for epsilon steps (see paper)
	 * Note: dt should be by definition dt = th - t0.
	 * Note: This function is for optimization purposes. I could assume that the whole Curve goes into computation, but then
	 *       the computing cost would grow, as Curve grows in time. Also, if some coefficients are not used
	 *       in rhs of the equation, we would carry out the unnecesary computations
	 *       (i.e. Jac Phi would be 0 and we would do 0 * C * r0).
	 *       There is one problem with this however: in rigorous code we need to check if the section is not crossed
	 *       between steps in Poincare map. Please see comments in JetSection / DDEPoincareMap classes.
	 *
	 * Param: out_u 			a finite-dimensional collection of data that is used to compute the map
	 * Param: admissible_order	will tell how high order of expansion is able to produce with data in u
	 */
	virtual void collectComputationData(
				const TimePointType& t0, const TimePointType& th, 	 	// input
				const RealType& dt, const CurveType& x, 	 			// input
				VariableStorageType& out_u, 							// output
				size_type& out_admissible_order) const = 0;				// output

	/**
	 * this is one of two computeDDECoefficients functions that needs to be implemented
	 * It takes as an input the proper set of variables representing input data to this map
	 * at time t0 (e.g. (appropriate subset of) output of collectComputationData())
	 * and makes use of it to compute recursively the Taylor expansion of the solution
	 * (coeffs) at t0. It also computes $\frac{\partial coeffs[i]}{\partial u}$ (in Du).
	 * It assumes that u contain appropriate data in u for this
	 * specific point t0, the user must assure this is true.
	 *
	 * THIS FUNCTION NEEDS TO BE IMPLEMENTED IN DERIVED CLASSES.
	 *
	 * Implementation note: the size of the ValueStorageType coeffs container
	 * is considered as a desired order. The user of the function should check
	 * if this is good size. The param out_admissible_order from collectComputationData()
	 * is used for this purpose.
	 */
	virtual void computeDDECoefficients(
				const RealType& t0, const ValueStorageType& u, 				// input
				ValueStorageType& coeffs, JacobianStorageType& Du) const = 0;	// output

	/**
	 * this is one of two computeDDECoefficients functions that needs to be implemented
	 * It takes as an input the proper set of variables representing input data to this map
	 * at time t0 (e.g. (appropriate subset of) output of collectComputationData())
	 * and makes use of it to compute recursively the Taylor expansion of the solution
	 * (coeffs) at t0. It assumes that u contain appropriate data in u for this
	 * specific point t0, the user must assure this is true. For non-developer users
	 * of the library it is possible to use versions of the function for specific CurveType.
	 *
	 * THIS FUNCTION NEEDS TO BE IMPLEMENTED IN DERIVED CLASSES.
	 *
	 * Implementation note: the size of the ValueStorageType coeffs container
	 * is considered as a desired order. The user of the function should check
	 * if this is good size. The param out_admissible_order from collectComputationData()
	 * is used for this purpose.
	 *
	 */
	virtual void computeDDECoefficients(
				const RealType& t0, const ValueStorageType& u, 	// input
				ValueStorageType& coeffs) const = 0;				// output

	/**
	 * computes recursively the Jet at time t for a given curve for a DDE of the form:
	 *
	 * 	x'(t) = F(t, x)
	 *
	 * Implementation note: the size of the ValueStorageType coeffs container
	 * is considered as a desired order. If order is 0 or order is too high for
	 * a given curve then the function should alter the order to highest possible.
	 * (use coeffs.resize(order) for this purpose)
	 *
	 * THERE IS DEFAULT IMPLEMENTATION OF THIS WITH collectComputationData() and
	 * computeDDECoefficients(..., u, ...) call to pure virtual function.
	 */
	virtual void computeDDECoefficients(
				const TimePointType& t0, const CurveType& x,
				ValueStorageType& coeffs) const {
		checkCurveDimension(x, "computeDDECoefficients");
		VariableStorageType u; size_type order;
		collectComputationData(t0, t0, RealType(0.), x, u, order);
		if (coeffs.size() == 0 || order + 1 < coeffs.size())
			coeffs.resize(order + 1);
		ValueStorageType v; convert(u, v);
		computeDDECoefficients(t0, v, coeffs);
	}

	/**
	 * Same as computeDDECoefficients(const TimePointType&, const CurveType&, JetType&),
	 * but also computes the partial derivative:
	 *
	 * 		Du[k] = \frac{\partial coeffs[k]}{\partial u}
	 *
	 * where u are all the variables used to evaluate the map.
	 * In u[j] is a d-dimensional set used in computation.
	 * In Du[k][j] is a matrix of dimension M(d,d) and of course it is
	 *
	 * 		Du[k][j] = \frac{\partial coeffs[k]}{\partial u[j]}
	 *
	 * It is done this way to allow Solvers to use this structure to reduce
	 * wrapping effect of interval arithmetics with the help of Lohner algorithm.
	 */
	virtual void computeDDECoefficients(
			const TimePointType& t0, const CurveType& x,
			ValueStorageType& coeffs, VariableStorageType& u, JacobianStorageType& Du) const {
		checkCurveDimension(x, "computeDDECoefficients");
		size_type order;
		collectComputationData(t0, t0, RealType(0.), x, u, order);
		if (coeffs.size() == 0 || order + 1 < coeffs.size())
			coeffs.resize(order + 1);
		ValueStorageType v;
		convert(u, v);
		computeDDECoefficients(t0, v, coeffs, Du);
	}

	/** this is for a current time in the solution. It provides basic implementation by call to the other function. */
	virtual VectorType operator()(const CurveType& x) const { return this->operator()(x.rightDomain(), x); }

	/** this is for a current time in the solution. It provides basic implementation by call to the other function. */
	virtual void computeDDECoefficients(const CurveType& curve, ValueStorageType& coeffs) const {
		computeDDECoefficients(curve.rightDomain(), curve, coeffs);
	}

	/** this is for a current time in the solution. It provides basic implementation by call to the other function. */
	virtual void computeDDECoefficients (
			const CurveType& curve,
			ValueStorageType& coeffs, VariableStorageType& u, JacobianStorageType& Du) const {
		computeDDECoefficients(curve.rightDomain(), curve, coeffs, u, Du);
	}

	/** virtual destructor for warning suppresion */
	virtual ~DDEBasicFunctionalMap() {}
protected:
	/** helper function to check. TODO: (NOT URGENT, FUTURE) make a switch to turn this off for speed.  */
	virtual void checkCurveDimension(CurveType const & x, std::string extra = "") const {
		if (x.dimension() != imageDimension()){
			std::ostringstream info;
			info << "FunctionalMap::" << extra << ": curve has incompatible dimension, ";
			info << "was " << x.dimension() << ", expecting " << imageDimension() << ".";
			throw std::domain_error(info.str());
		}
	}
};




/**
 * TODO: rewrite DOCS
 * this is the interface of a Functional map, that can be evaluated
 * only by "querying" the curve about its jets at various points in the past.
 * it can produce jet at a given time t together with partial derivatives
 * w.r.t. to internal variables from curve used in the computation by
 * recurrent formula for a functional equation:
 *
 * 		x'(t) = f(t, x)
 *
 * where x is the solution curve in \R^d defined for times smaller than
 * t (up to some maximal past time usually).
 *
 * A good example of such a functional map is discrete delay case,
 * e.g. F(t, x) = f(t, x(t), x(t-tau)), f : \R \times \R^{2d}->\R
 * (easily extended to more delays). (Probably) the functionals
 * that include integral over the past: F(t, x) = f(t, \int_{t-tau}^{t} x(s) ds)
 * will also fall into this category.
 *
 * The concrete implementations will be available later.
 *
 * Implementation note: when implementing concrete class of
 * this interface, do not forget to put
 * 	 using BaseClass::operator();
 *	 using BaseClass::computeDDECoefficients;
 * in the public part of your class to gain access to default
 * implementations of some functions.
 *
 * TODO: (NOT URGENT, FUTURE, RETHINK) it would be better not to pass containers to computeDDECoefficients, but maybe iterators to given types? It could be used to compute in a given place, i.e. at an already created Jet?
 */
template<typename SolutionCurveSpec, typename JetSpec=typename SolutionCurveSpec::JetType>
class DDERigorousFunctionalMap : public DDEBasicFunctionalMap<SolutionCurveSpec, JetSpec> {
public:
	typedef DDEBasicFunctionalMap<SolutionCurveSpec, JetSpec> BaseClass;
	typedef typename BaseClass::CurveType CurveType;
	typedef typename BaseClass::JetType JetType;
	typedef typename BaseClass::RealType RealType;
	typedef typename BaseClass::TimePointType TimePointType;
	typedef typename BaseClass::MatrixType MatrixType;
	typedef typename BaseClass::VectorType VectorType;
	typedef typename BaseClass::ScalarType ScalarType;
	typedef typename BaseClass::size_type size_type;
	typedef typename BaseClass::DataType DataType;
	typedef typename BaseClass::VariableStorageType VariableStorageType;
	typedef typename BaseClass::ValueStorageType ValueStorageType;
	typedef typename BaseClass::JacobianStorageType JacobianStorageType;

	using BaseClass::convert;
	using BaseClass::imageDimension;
	using BaseClass::operator();
	using BaseClass::collectComputationData;
	using BaseClass::computeDDECoefficients;

	/**
	 * Extra method needed by the rigorous code
	 *
	 * Collects all data from curve necessary for the computation of the value.
	 * We use  collection of DataType items as the output, we ASSUME that SetType elements are
	 * used to describe function by CurveType. Therefore, in the output, we could
	 * relate how inputs (VariableStorage) corresponds to outputs (ValueStorage)
	 * mainly by Jacobian of the map at Variables (JacobianStorageType).
	 *
	 * Then, the main function that needs to be implemented in the Map
	 * are those depending on the VariableStorageType. Other can have
	 * default implementations and are for users convenience.
	 *
	 * THIS FUNCTION NEEDS TO BE IMPLEMENTED IN DERIVED CLASSES.
	 *
	 * Note: dt might seem redundant here, but for some reasons it might be helpful, i.e. for epsilon steps (see paper)
	 * Note: dt should be by definition dt = th - t0.
	 * Note: This function is for optimization purposes. I could assume that the whole Curve goes into computation, but then
	 *       the computing cost would grow, as Curve grows in time. Also, if some coefficients are not used
	 *       in rhs of the equation, we would carry out the unnecesary computations
	 *       (i.e. Jac Phi would be 0 and we would do 0 * C * r0).
	 *       There is one problem with this however: in rigorous code we need to check if the section is not crossed
	 *       between steps in Poincare map. Please see comments in JetSection / DDEPoincareMap classes.
	 *
	 * Param: out_u 			a finite-dimensional collection of data that is used to compute the map
	 * Param: out_encl 			first entries corresponds to enclosures over [t0, t0+dt] of the variables in u,
	 * 							then, it might contain more entries, relevant to computation. For example higher
	 * 							order enclosures used in Taylor solver.
	 * Param: admissible_order	will tell how high order of expansion is able to produce with data in u
	 */
	virtual void collectComputationData(
				const TimePointType& t0, const TimePointType& th, 	 				// input
				const RealType& dt, const CurveType& x, 	 						// input
				VariableStorageType& out_u, ValueStorageType& out_encl,				// output
				size_type& out_admissible_order) const = 0;							// output

	/** virtual destructor for warning suppresion */
	virtual ~DDERigorousFunctionalMap() {}
protected:
	using BaseClass::checkCurveDimension;
};


} // namespace ddes
} // namespace capd


#endif /* _CAPD_DDES_FUNCTIONALMAP_H_ */
