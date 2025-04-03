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

#ifndef _CAPD_DDES_DDESSLIB_H_
#define _CAPD_DDES_DDESSLIB_H_

// those are the common things both non-rigorous and rigorous code uses
#include <capd/ddes/DDECommon.hpp>
#include <capd/ddes/storage/storage.h>
#include <capd/ddes/FunctionalMap.hpp>
#include <capd/ddes/DDEJetSection.hpp>
#include <capd/ddes/SampleEqns.hpp>

// those are for nonrigorous code
#include <capd/ddes/DDEPiecewisePolynomialCurve.hpp>
#include <capd/ddes/BasicDiscreteDelaysFunctionalMap.hpp>
#include <capd/ddes/DDENonrigorousTaylorSolver.hpp>
#include <capd/ddes/DDEBasicPoincareMap.hpp> // TODO: (!!!URGENT) PRAWDOPODOBNIE USUNAC JAK SKONCZE REFACTOR! (PoincareMap bedzie sciagac dane z tego)
#include <capd/ddes/DDEBasicTimeMap.hpp>

// those are for rigorous code.
#include <capd/ddes/DDEForwardTaylorCurvePiece.hpp>
#include <capd/ddes/DDESolutionCurve.hpp>
#include <capd/ddes/DiscreteDelaysFunctionalMap.hpp>
#include <capd/ddes/DDETaylorSolver.hpp>
#include <capd/ddes/DDEPoincareMap.hpp>

namespace capd {

namespace ddes {


// TODO: (IMPORTANT) DOCS!
// TODO: (IMPORTANT) TESTS!
template<
	typename EqSpec,
	typename MatrixSpec = capd::DMatrix,
	typename VectorSpec = capd::DVector
>
class NonrigorousSetup{
public:
	typedef EqSpec Eq;
	//typedef typename Eq::ParamType ParamType; // TODO: removed from here, as in the basic setup user should control them themself
	typedef VectorSpec Vector;
	typedef MatrixSpec Matrix;
	typedef typename Matrix::ScalarType Scalar;
	typedef typename Matrix::ScalarType Real; // TODO: (FUTURE) Rethink? What if scalar is Complex?
	//typedef typename Eq::ParamsVectorType ParamsVector; // TODO: removed from here, as in the basic setup user should control them themself
	typedef capd::ddes::DiscreteTimeGrid<Real> Grid;  ///< important class, defining the computation grid to produce $t_i$ points
	typedef typename Grid::TimePointType TimePoint;	  ///< important class, defining the grid points $t_i = ih$
	typedef capd::ddes::GenericJet<TimePoint, Vector, Vector, Matrix> Jet;   ///< this is to be used in C1Solution, you will probably not use this directly
	typedef capd::ddes::GenericJet<TimePoint, capd::ddes::VectorWithJacData<Vector, Matrix>, Vector, Matrix> C1Jet;  ///< this is to be used in C1Solution, you will probably not use this directly
	typedef capd::ddes::DDEPiecewisePolynomialCurve<Grid, Jet> Solution;     ///< basic class to store solutions
	typedef capd::ddes::DDEPiecewisePolynomialCurve<Grid, C1Jet> C1Solution; ///< a class to store solutions when doing C^1 computations (e.g. computing Jacobian of a Poincare Map)
	typedef Solution Segment; 		///< alias for Solution
	typedef C1Solution C1Segment; 	///< alias for C1Solution
	typedef Jet CurvePiece;
	typedef capd::ddes::BasicDiscreteDelaysFunctionalMap<Eq, C1Solution> C1DDEq;
	typedef capd::ddes::BasicDiscreteDelaysFunctionalMap<Eq, Solution> DDEq;
	typedef capd::ddes::DDENonrigorousTaylorSolver<C1DDEq> C1Solver;
	typedef capd::ddes::DDENonrigorousTaylorSolver<DDEq> Solver;
	typedef typename C1Solver::VariableStorageType Variables;
	typedef typename C1Solver::JacobianStorageType Jacobians;
	typedef typename C1Solver::ValueStorageType Values;
	typedef typename C1Solver::size_type size_type;
	typedef capd::ddes::DDEJetSection<C1Solution> C1Section;
	typedef typename C1Section::JetType C1SecJet;
	typedef capd::ddes::DDEBasicPoincareMap<C1Solver, C1Section> C1PoincareMap;
	typedef capd::ddes::DDEJetSection<Solution> Section;
	typedef capd::ddes::DDEBasicPoincareMap<Solver, Section> PoincareMap;

	typedef capd::DMap CAPDMap;
};


// TODO: (IMPORTANT) DOCS!
// TODO: (IMPORTANT) TESTS!
template<
	typename EqSpec,
	typename MatrixSpec = capd::IMatrix,
	typename VectorSpec = capd::IVector,
	typename PoliciesSpec=capd::dynset::C11Rect2Policies
>
class RigorousSetup{
public:
	typedef EqSpec Eq;
	// typedef typename Eq::ParamType ParamType; // TODO: removed from here, as in the basic setup user should control them themself
	typedef VectorSpec Vector;
	typedef MatrixSpec Matrix;
	typedef typename Matrix::ScalarType Scalar;
	typedef typename Matrix::ScalarType Real; // TODO: (FUTURE) Rethink? What if scalar is Complex?
	// typedef typename Eq::ParamsVectorType ParamsVector; // TODO: removed from here, as in the basic setup user should control them themself
	typedef PoliciesSpec Policies;
	typedef capd::ddes::SharedDoubleton<Matrix, Policies> SetType;
	typedef capd::ddes::DDESolutionCurve<SetType> Solution;
	typedef typename Solution::GridType Grid;
	typedef typename Solution::TimePointType TimePoint;
	typedef typename Solution::CurvePieceType CurvePiece;
	typedef capd::ddes::DiscreteDelaysFunctionalMap<Eq, Solution> DDEq;
	typedef capd::ddes::DDETaylorSolver<DDEq> Solver;
	typedef typename Solver::VariableStorageType Variables;
	typedef typename Solver::JacobianStorageType Jacobians;
	typedef typename Solver::ValueStorageType Values;
	typedef typename Solver::size_type size_type;
	typedef capd::ddes::DDEJetSection<Solution> Section;
	typedef typename Section::JetType SectionJet;
	typedef capd::ddes::DDEPoincareMap<Solver, Section> PoincareMap;

	typedef capd::IMap CAPDMap;
};

} // namespace ddes;

} // namespace capd;

#endif /* _CAPD_DDES_DDESSLIB_H_ */
