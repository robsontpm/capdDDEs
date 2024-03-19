/*
 * Common definitions. As DDE code library is not yet standalone
 * and we do not have the clever parser of the RHS as in CAPD codes
 * we need to redefine some types and functions.
 *
 * Also, we define some helper functions for the output
 */
#ifndef EXAMPLES_ROSSLER_ODE_VS_DDE_CODE_SETUP_H_
#define EXAMPLES_ROSSLER_ODE_VS_DDE_CODE_SETUP_H_

#include <iostream>
#include <iomanip>
#include <capd/capdlib.h>
#include <capd/ddes/ddeslib.h>
#include <capd/ddeshelper/ddeshelperlib.h>
#include <capd/ddeshelper/DDEHelperRigorous.hpp>
#include "equation.h"

using namespace std;
using namespace capd;

/**
 * we will use Rossler DDE with interval parameters.
 * This is the definition of f from x'(t) = f(x(t), x(t-tau)).
 */
typedef Rossler<interval> Eq;

/** this contains all necessary ingredients to do rigorous numerics in DDEs */
typedef capd::ddeshelper::RigorousHelper<Eq, 1, IMatrix, IVector> Setup;

// below are just renaming for shorter class names
typedef Setup::Grid Grid;				// grid for the solutions in the C^n_p space.
typedef Setup::DDEq DDEq;				// this is F from the abstract formulation x'(t) = F(x_t), that contains the information on delays. In our case F(x_t) := f(x_t(0), x_t(-tau)).
typedef Setup::Solution DDESolution;	// basic set type for DDEs
typedef Setup::Solver DDESolver;		// semidynamical system
typedef Setup::PoincareMap DDEPoincare;	// Poincare map for the semidynamical system
typedef Setup::Section DDESection;		// basic section for DDEs

#endif /* EXAMPLES_ROSSLER_ODE_VS_DDE_CODE_SETUP_H_ */
