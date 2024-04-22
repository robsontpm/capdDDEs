/*
 * Here we put only some constants that are common for both
 * non-rigorous and rigorous phases of the computations.
 */
#ifndef EXAMPLES_MACKEY_GLASS_STABLE_PERIODIC_COMMON_PARAMETRS_H_
#define EXAMPLES_MACKEY_GLASS_STABLE_PERIODIC_COMMON_PARAMETRS_H_

namespace PARAMS {
	// params of the (p,n)-representation
	// they MUST be the same for non-rigorous and rigorous computations
	// as non-rigorous must produce data readable by rigorous code
	// the meaning of the parameters are described in our papers, see references elsewhere.
	const int p = 128;
	const int n = 4;

	// setup parameters of the Mackey-Glass equation
	// they will be used in both codes as default,
	// but you might change those in the rig-setup.h to their interval representations if needed.
	// in case those parameters are used for both codes, then they should be representable real numbers
	// otherwise the proof will be not for those values but for representable values close to those.
	// Names of the parameters are as in the article: http://www.scholarpedia.org/article/Mackey-Glass_equation
	const double TAU { 2.0 };
	const double BETA { 2.0 };
	const double GAMMA { 1.0 };
	const double N { 6.0 }; // you might also try n=8, as in FOCM 2018 paper.

	// setup for the nonrigorous procedure to find the candidate
	// the method is as follows:
	// 1. integrate constant solution x_0(t) = v0 for all t \in [-tau, 0] for NUM_ITERATES full delays
	// 2. approach the selected section from the last segment.
	// 3. for the segment on the section, compute Poincare map NUM_POINCARE times
	// 4. apply Newton method on the poincare map for the candidate for at least NUM_NEWTON steps, or until difference ||x0 - P(x0)||_sup < EPS_NEWTON
	const double v0 = 1.1;			 // value of the initial constant solution segment to start nonrigorous computations
	const int NUM_ITERATES =  10;	 // number of full delays to integrate the constant solution, before trying to do poincare map
	const int NUM_POINCARE =  10;	 // number of poincare images to compute before proceeding to
	const int NUM_NEWTON =    10;	 // number of Newton method steps to apply to find 0 of P(x0) - x0.
	const double EPS_NEWTON = 6e-12; // or stop if ||x0 - P(x0)||_sup < EPS_NEWTON. Rule of thumb: it should be greater than machine_precision * dimension M() of the vector, so usually around 10^-13 or at most 10.

	const int REQUIRED_STEPS = -(n+1);	 // number of full steps to be made before attempting to go back to section. Negative value means full delays. This should be made ,,large enough'' (i.e. greater than n+1), but you might need greater values if you are tracing really long solution.
}

#endif /* EXAMPLES_MACKEY_GLASS_STABLE_PERIODIC_COMMON_PARAMETRS_H_ */
