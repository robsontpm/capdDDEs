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

#ifndef _CAPD_DDES_BASICDISCRETEDELAYSFUNCTIONALMAP_H_
#define _CAPD_DDES_BASICDISCRETEDELAYSFUNCTIONALMAP_H_

#include <capd/ddes/DDECommon.h>
#include <capd/ddes/FunctionalMap.h>
#include <capd/basicalg/TypeTraits.h>
#include <stdexcept>
#include <string>
#include <capd/fadbad/fadbad.h>
#include <capd/fadbad/tadiff.h>
#include <capd/fadbad/fadiff.h>
// DEV NOTE: here we do not need #include <capd/dynsys/filibfadbad.h>, as it is dor
// DEV NOTE: interval code, and BasicFunctionalMap is just for no nrigorous computations

namespace capd{
namespace ddes{

/**
 * A class to represent right hand side of the DDE of the form:
 *
 *     $x' = f(t, x(t), x(t-tau_1), ..., x(t-tau_m))$
 *
 * NOTE: This is a demo version of what I want to achieve with the DDEs code interface
 *       I assume it should be similar to ODEs in original CAPD.
 *       I assume that any RHS map f of the equation $x'(t) = f(x_\tau)$
 *       should be representable with similar interface (i.e. the function
 *       accepts some SolutionCurve and can produce Jet at current time in the
 *       solution defined with the equation (function computeDDECoefficients)
 *       I do not assume anything about SolutionCurve except, it can produce
 *       Forward Taylor Jets of itself at a given point (that is used in the equation)
 *       and that it can return its value at a current time.
 *
 *       In the implementation below, I assume that I can ask about SolutionCurve and its Jet
 *       at $t-\tau_i$, where $t$ is current time. Then I construct an AD jet of
 *       $x(t-\tau_i)$ for all necessary $i$'s and pass sequentially through
 *       the procedure similar to that of CAPD computeODECoefficients().
 *
 * NOTE: I do not know if integral equations (i.e. $x'(t) = F(\int_{t-\tau}^t g(x(s)) ds$)
 * 		 could be described in this manner but I hope it could be done, at least if g is
 * 		 linear in x (e.g. g = Id).
 * NOTE: in fact, functional map could have operators that are templates, but we assume that
 *       the argument curve might have some requirements when computing, for example
 *       it could determine only the delays at the specific grid points (because of the loss of regularity,
 *       we are not able to always give jets over some general intervals, see new notes).
 *       Therefore, we supply the Map with appropriate CurveType that it can manipulate.
 * NOTE: This version is used only in non-rigorous computations!
 * NOTE: The @see capd::ddes::NonrigorousSetup and @see capd::ddeshelper::NonrigorousHelper classes
 *       defines a concrete version of this class for you, so you don't need to worry
 *       about the template parameters and such.
 *
 * Template parameters:
 * @param FinDimMapSpec this is a specification of f in r.h.s., we assume $f : \R^{d + m*d} \to \R^d$
 * @param SolutionCurveSpec the representation of the solution curve, the most important thing that this class depends on
 * @param JetSpec the default value should always be ok
 */
template<
	typename FinDimMapSpec,
	typename SolutionCurveSpec,
	typename JetSpec=typename SolutionCurveSpec::JetType
>
class BasicDiscreteDelaysFunctionalMap : public DDEBasicFunctionalMap<SolutionCurveSpec, JetSpec> {
public:
	typedef DDEBasicFunctionalMap<SolutionCurveSpec, JetSpec> BaseClass; 	///< User Friendly Renaming. It is very convenient e.g. when calling the constructor of base class!
	typedef BasicDiscreteDelaysFunctionalMap<FinDimMapSpec, SolutionCurveSpec, JetSpec> Class; ///< User Friendly Renaming. It is very convenient when declaring functions working on this class or returning this class. It's like Java and I like it :)
	typedef SolutionCurveSpec CurveType;		///< User Friendly Renaming. Data type that can be used to evaluate the functional map on
	typedef JetSpec JetType;					///< User Friendly Renaming. Type for the Jet of the solution
	typedef FinDimMapSpec MapType;				///< User Friendly Renaming. Type for the Jet of the solution
	typedef typename BaseClass::RealType RealType;			///< User Friendly Renaming. In theory Scalar might be something like Complex, so this might differ from ScalarTtpe, but for now we support only ScalarType here.
	typedef typename BaseClass::TimePointType TimePointType;///< User Friendly Renaming. Time point type for the grid used in computations. See papers for the idea of the time grid.
	typedef typename BaseClass::MatrixType MatrixType;		///< User Friendly Renaming.
	typedef typename BaseClass::VectorType VectorType;		///< User Friendly Renaming.
	typedef typename BaseClass::ScalarType ScalarType;		///< User Friendly Renaming. Scalar type for the Matrix and Vector types.
	typedef typename BaseClass::size_type size_type;		///< User Friendly Renaming. Usually some unsigned integer.
	typedef typename BaseClass::DataType DataType;			///< User Friendly Renaming. DataType to be used to represent the chunks of the finite representation of Curve. Basically it is VectorType, but can differ in some instances. In rigorous code it is more important.
	typedef typename BaseClass::ValueStorageType ValueStorageType;			///< User Friendly Renaming. The data type used by this class to store Raw data needed for computation. This will usually be a collection (e.g. std::vector) of VectorType
	typedef typename BaseClass::VariableStorageType VariableStorageType;	///< User Friendly Renaming. The data type used by this class to store Raw data needed for computation. This will usually be a collection (e.g. std::vector) of DataType
	typedef typename BaseClass::JacobianStorageType JacobianStorageType;	///< User Friendly Renaming. The data type used by this class to store Raw data for C^11 computation.  This will usually be a collection (e.g. std::vector) of MatrixType
	/** storage type for discrete number of delays */
	typedef std::vector< TimePointType > DelayStorageType;
	/** iterator over this map will iterate over delays! */
	typedef typename DelayStorageType::const_iterator const_iterator;
	/** iterator over this map will iterate over delays! */
	typedef typename DelayStorageType::iterator iterator;
	/** standard rebind of the types as in many C++ std classes */
    template <typename OtherSolutionSpec, typename OtherJetSpec=typename OtherSolutionSpec::JetType>
    struct rebind { typedef BasicDiscreteDelaysFunctionalMap<FinDimMapSpec, OtherSolutionSpec, OtherJetSpec> other; };

    /** conversion between e.g. rigorous and non-rigorous version */
	template<typename OtherDiscreteDelaysFunctionalMap>
	BasicDiscreteDelaysFunctionalMap(const OtherDiscreteDelaysFunctionalMap& orig): m_map(orig.m_map), m_delays(orig.m_delays){}
	/** standard copy contructor */
	BasicDiscreteDelaysFunctionalMap(const BasicDiscreteDelaysFunctionalMap& orig): m_map(orig.m_map), m_delays(orig.m_delays){}
	/**
	 * this is the basic version, you should supply the finite map and the delays.
	 *
	 * You might do it this way (you need to use your own data types!), using initializers { }:
	 *
	 * 		MyMap f;
	 * 		DiscreteDelaysFunctionalMap<MyMap, SolutionCurve>(f, {tau1, tau2});
	 *
	 * Please note that the constructors check if you supplied enough delays to match the signature of the MapType!
	 * @see SampleEqns.h for some examples how to write your own MapType!
	 *
	 * If the number of delays do not match the MapType, then an std::error will be thrown.
	 */
	BasicDiscreteDelaysFunctionalMap(const MapType& f, std::vector<TimePointType> delays): m_map(f), m_delays(delays){ this->checkMapSignature(); }
	/**
	 * this uses default constructor of the MapType and assume no delays!
	 * @see BasicDiscreteDelaysFunctionalMap(const MapType& f, std::vector<TimePointType> delays)
	 */
	BasicDiscreteDelaysFunctionalMap(){ this->checkMapSignature(); }
	/**
	 * this uses default constructor of the MapType and assume there is just one delay!
	 * @see BasicDiscreteDelaysFunctionalMap(const MapType& f, std::vector<TimePointType> delays)
	 */
	BasicDiscreteDelaysFunctionalMap(const TimePointType& tau) { m_delays.push_back(tau); this->checkMapSignature(); }
	/**
	 * this uses the given map f and assume there is no delay!
	 * @see BasicDiscreteDelaysFunctionalMap(const MapType& f, std::vector<TimePointType> delays)
	 */
	BasicDiscreteDelaysFunctionalMap(const MapType& f): m_map(f) { this->checkMapSignature(); }
	/**
	 * this uses the given map f and assume there is just one delay!
	 * This is for convenience for the most users.
	 * @see BasicDiscreteDelaysFunctionalMap(const MapType& f, std::vector<TimePointType> delays)
	 */
	BasicDiscreteDelaysFunctionalMap(const MapType& f, const TimePointType& tau): m_map(f) { m_delays.push_back(tau); this->checkMapSignature(); }
	/**
	 * this is a constructor in an old C like fashion. Maybe it should be deprecated?
	 * n is the size of delays.
	 * @see BasicDiscreteDelaysFunctionalMap(const MapType& f, std::vector<TimePointType> delays)
	 */
	BasicDiscreteDelaysFunctionalMap(const MapType& f, int n, TimePointType delays[]): m_map(f){ for (size_type i = 0; i < n; ++i) m_delays.push_back(delays[i]); this->checkMapSignature(); }
	/** standard */
	virtual ~BasicDiscreteDelaysFunctionalMap() {}

	/** iterator over delays */
	iterator begin() { return m_delays.begin(); }
	/** iterator over delays */
	iterator end() { return m_delays.end(); }
	/** iterator over delays */
	const_iterator begin() const { return m_delays.begin(); }
	/** iterator over delays */
	const_iterator end() const { return m_delays.end(); }

	/** output dimension of the internal map */
	size_type imageDimension() const { return m_map.imageDimension(); }
	/** input dimension of the internal map, x(t) is one, x(t-tau_1) is another, etc.. (thus the formula) */
	size_type dimension() const { return (m_delays.size() + 1) * m_map.imageDimension(); }
	/** number of delays */
	size_type delaysCount() const { return m_delays.size(); }

	/**
	 * Implementation of the evaluation of this function on the given curve.
	 * @see BaseClass documentation for more details.
	 */
	VectorType operator()(const TimePointType& t, const CurveType& x) const;

	/**
	 * This is an implementation of a virtual func. that is needed
	 * by Taylor-type integrators. See Interface docs for information
	 * on what is expected to be returned in given variables.
	 *
	 * Input is self explanatory: t0 = current time, th = t0 + h = next time,
	 * dt = h = the step size, x - the solution curve.
	 * Those are kind of redundant, but they might be useful in some other scenarios,
	 * like computing epsilon steps, unbounded delays, etc.
	 * The curve x is not necessary the solution segment x_t, it might be the
	 * whole solution on a big interval, but it should span at least [t0 - max delay, t0].
	 *
	 * The output (t = t0 here):
	 * out_u = $(x(t), x(t-\tau_1), .., x(t-\tau_m), x^[1](t-\tau_1), ..., x^[1](t-\tau_m), x^[2](t-\tau_1), ..., x^[k](t-\tau_1), ... )$
	 * out_admissible_order = the maximal order of a jet you are supposed to produce in your computations
	 * you can use jests to out_admissible_order - 1, i.e. the order k in the definition of u
	 * above goes up to out_admissible_order - 1.
	 */
	virtual void collectComputationData(
					const TimePointType& t0, const TimePointType& th,
					const RealType& dt, const CurveType& x,
					VariableStorageType& out_u,
					size_type& out_admissible_order) const;
	/**
	 * Implementation of virtual func. See Interface docs.
	 * Basically it computes Taylor coefficients of the solution at t0
	 */
	virtual void computeDDECoefficients(
					const RealType& t0, const ValueStorageType& u,
					ValueStorageType& coeffs) const;

	/**
	 * see Interface docs. This is more involved version, that also
	 * computes derivative of the step map with respect of initial data.
	 * This might be used by Lohner algorithm in the set.move() to select
	 * good coordinate frame for the next set.
	 *
	 * The steps are basically the same as in
	 * @see virtual void computeDDECoefficients(const RealType&, const ValueStorageType&, ValueStorageType&),
	 * but we do some extra steps. The comments from that procedure apply to this.
	 */
	virtual void computeDDECoefficients(
				const RealType& t0, const ValueStorageType& u,
				ValueStorageType& coeffs, JacobianStorageType& Du) const;

	/**
	 * returns maximal delay.
	 * Note: running time O(m), m - number of delays.
	 */
	TimePointType getMaxDelay();
	/**
	 * Returns a segment that spans from  to max delay.
	 * The segment represents the constant solution equal to the given vector over whole interval.
	 * Useful to start computation with something.
	 */
	CurveType makeCompatibleSegment(VectorType v={});

	using BaseClass::operator();
	using BaseClass::computeDDECoefficients;
protected:
	MapType m_map;
	DelayStorageType m_delays;

	/**
	 * For internal sanity checks.
	 * checks if the given dimension is compatible with the output of this map.
	 * Throws std::logic_error with a verbose message if not.
	 */
	void checkDimension(size_type d) const;
	/**
	 * For internal sanity checks.
	 * Checks if the finite dimensional map is compatible with all the delays given.
	 * Throws std::logic_error with a verbose message if not.
	 */
	void checkMapSignature() const;

};


} // namespace ddes
} // namespace capd


#endif /* _CAPD_DDES_BASICDISCRETEDELAYSFUNCTIONALMAP_H_ */
