/*
 * Common definitions. As DDE code library is not yet standalone
 * and we do not have the clever parser of the RHS as in CAPD codes
 * we need to redefine some types and functions.
 *
 * Also, we define some helper functions for the output
 */
#ifndef EXAMPLES_MACKEY_GLASS_STABLE_PERIODIC_NONRIG_SETUP_H_
#define EXAMPLES_MACKEY_GLASS_STABLE_PERIODIC_NONRIG_SETUP_H_

#include <iostream>
#include <capd/capdlib.h>
#include <capd/ddes/ddeslib.h>
#include <capd/ddeshelper/ddeshelperlib.h>
#include <capd/ddeshelper/DDEHelperNonrigorous.hpp>
#include "common-params.h"

using namespace std;

// we say that the system computes on Intervals and has Intervals as parameters
typedef capd::ddes::MackeyGlass<double, double> Eq;

// the number one (1) is the number of delays
typedef capd::ddeshelper::NonrigorousHelper<Eq, 1, capd::DMatrix, capd::DVector> DDEs;

// setup parameters of the system
// by default, we take the values of the double parameters
// defined for both systems in common-parametrs.h, but they
// might be changed here, for example to rigorous bounds on the parameters
// if they are non-representable real numbers.
namespace PARAMS {
	// last one (i.e. 1.0) is the delay used in computations as the length of the delay
	// we scale the parameters of the equations by PARAM_TAU. This is a standard thing.
	const typename DDEs::ParamsVector params {GAMMA * TAU, BETA * TAU, N, 1.0};
}

#endif /* EXAMPLES_MACKEY_GLASS_STABLE_PERIODIC_NONRIG_SETUP_H_ */
