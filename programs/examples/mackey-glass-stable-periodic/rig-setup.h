/*
 * Common definitions. As DDE code library is not yet standalone
 * and we do not have the clever parser of the RHS as in CAPD codes
 * we need to redefine some types and functions.
 *
 * Also, we define some helper functions for the output
 */
#ifndef EXAMPLES_MACKEY_GLASS_STABLE_PERIODIC_RIG_SETUP_H_
#define EXAMPLES_MACKEY_GLASS_STABLE_PERIODIC_RIG_SETUP_H_

#include <iostream>
#include <capd/capdlib.h>
#include <capd/ddes/ddeslib.h>
#include <capd/ddeshelper/ddeshelperlib.h>
#include <capd/ddeshelper/DDEHelperRigorous.hpp>
#include "common-params.h"

using namespace std;

// common setup for the DDE integrator machinery,
// we will use old capd intervals, as the autodiff code is not ready for new intervals
// especially, the ln() function used in the definition of the MackeyGlass in SampleEqns.h
typedef capd::intervals::Interval<double, capd::rounding::DoubleRounding>  Interval;

// we say that the system computes on Intervals and has Intervals as parameters
typedef capd::ddes::MackeyGlass<Interval, Interval> Eq;

// the number one (1) is the number of delays
typedef capd::ddeshelper::RigorousHelper<Eq, 1> RigSetup;

// this helps output nice comparisons for high-dimensional sets
typedef capd::ddeshelper::DDECompareHelper<RigSetup::Vector> Comparator;

// setup parameters of the system
// by default, we take the values of the double parameters
// defined for both systems in common-parametrs.h, but they
// might be changed here, for example to rigorous bounds on the parameters
// if they are non-representable real numbers.
namespace RIGPARAMS {
	const Interval TAU { PARAMS::TAU };
	const Interval BETA { PARAMS::BETA };
	const Interval GAMMA { PARAMS::GAMMA };
	const Interval N { PARAMS::N };

	// last one (i.e. 1.0) is the delay used in computations as the length of the delay
	// we scale the parameters of the equations by PARAM_TAU. This is a standard thing.
	const RigSetup::ParamsVector params {GAMMA * TAU, BETA * TAU, N, 1.0};
}

#endif /* EXAMPLES_MACKEY_GLASS_STABLE_PERIODIC_RIG_SETUP_H_ */
