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

	// setup for the nonrigorous procedure to find the candidate
	// the method is as follows:
	// 1. integrate constant solution x_0(t) = v0 for all t \in [-tau, 0] for NUM_ITERATES full delays
	// 2. approach the selected section from the last segment.
	// 3. for the segment on the section, compute Poincare map NUM_POINCARE times
	// 4. apply Newton method on the poincare map for the candidate for at least NUM_NEWTON steps, or until difference ||x0 - P(x0)||_sup < EPS_NEWTON
	const double v0 = 1.1;			 // value of the initial constant solution segment to start nonrigorous computations
	const int NUM_ITERATES =  10;	 // number of full delays to integrate the constant solution, before trying to do poincare map
	const int NUM_POINCARE =  10;	 // number of poincare images to compute before proceeding to
	const int NUM_NEWTON =     0;	 // number of Newton method steps to apply to find 0 of P(x0) - x0.
	const double EPS_NEWTON = 6e-12; // or stop if ||x0 - P(x0)||_sup < EPS_NEWTON. Rule of thumb: it should be greater than machine_precision * dimension M() of the vector, so usually around 10^-13 or at most 10.
}

#endif /* EXAMPLES_MACKEY_GLASS_STABLE_PERIODIC_NONRIG_SETUP_H_ */
