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

	// number of full steps to be made before attempting to go back to section.
	// Negative value means full delays. This should be made ,,large enough''
	// (i.e. greater than n+1), but you might need greater values if you are tracing really long solution.
	const int REQUIRED_STEPS = -(n+1);

	// the radius of the initial set in the dominant eigendirection.
	// Others will be guessed from that radius by heuristic procedure.
	const double INITIAL_RADIUS = 1e-4;
}

#endif /* EXAMPLES_MACKEY_GLASS_STABLE_PERIODIC_COMMON_PARAMETRS_H_ */
