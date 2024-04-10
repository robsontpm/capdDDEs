/* * TODO: Docs */// comment or uncomment the following for extra options// DDES_ALLOW_SYSTEM 		when uncommented, the program can produce some pictures from the proof,//							if you only interested in the proof of (2) then you can uncomment it// EXTRA_OUTPUT				uncommenting this will output more information to the standard output#define DDES_ALLOW_SYSTEM//#define EXTRA_OUTPUT// ====================== DO NOT MODIFY ANYTHING BELOW THIS LINE ===========================================#include "nonrig-setup.h"DDEs::Vector find_candidate(DDEs& setup, DDEs::Vector v0);
int main(int, char**){	DDEs setup(PARAMS::params, PARAMS::p, PARAMS::n, PARAMS::REQUIRED_STEPS);	DDEs::Vector inital = find_candidate(setup, { PARAMS::v0 });	return 0;}DDEs::Vector find_candidate(DDEs& setup, DDEs::Vector v0){	const DDEs::Grid& grid = setup.grid();			// for easier access to grid	auto plot_dt = double(grid(1))/10.; 			// used in plots.	int steps = PARAMS::NUM_ITERATES * setup.p(); 	// how many full steps to do before starting refinement process, we convert parameter to full intervals!	auto zero = grid(0);		// t=0, as a time point on the grid	auto tau = grid(setup.p()); // tau as a time point on the grid	cout << "phase 0 start - refine by iteration" << endl;	// make a constant x(s) = v0 solution segment on [-tau, 0] with an order n of jets set in the setup.	// We will use it to convert between Vector (no structure, just \R^M) and Solution (with structure, jets, grid, etc.)	// representation of various segments, when necessary.	DDEs::Solution candidate = setup.makeSegment(v0);	auto trajectory = setup.integrate(steps, candidate, candidate);	// NOTE: with the helper functions such as integrate, poincare, make***(), etc.	//       you can use auto keyword in C++ 11 and newer to stop caring about the typenames!	capd::ddeshelper::plot_value("phase-0-initial--", plot_dt, trajectory, false);	DDEs::Solution s(grid, -tau, zero, setup.n(), { 0.0 }); //you may use also setup.makeSegment(); instead. This is just to show how to make segments by hand.	s.value(zero) = { -1.0 };			// 0-th coefficient at t=0	s.value(-tau) = { +1.0 };			// 0-th coefficient at t=-tau	DDEs::JetSection section(s, 0.0); 	// this defines the section S = {x : x(0) == x(-tau)}.	trajectory = setup.poincare(section, candidate, candidate);	// we are re-using candidate as the output 3rd parameter, i.e. P(candidate)	// we try to refine the periodic orbit by interation (we hope it is attracting)	capd::ddeshelper::plot_value("phase-0-candidate-on-section--", plot_dt, candidate, false);	for (int i = 0; i < PARAMS::NUM_POINCARE; ++i)		setup.poincare(section, candidate, candidate);	// now, we test how good the candidate is up to now.	// We are after stable solution, so iteration should produce quite good results.	DDEs::Solution segment = setup.makeSegment(); // make a const=0 function wth a proper structure.	trajectory = setup.poincare(section, candidate, segment);  // segment = P(candidate)	capd::vectalg::EuclLNorm<DDEs::Vector, DDEs::Matrix> euclNorm; // we use simple norm to assess the quality of the candidate	// for now, the arithmetic operations are not implemented in Solution, even for the same structures.	// so we just convert to vector to compute standard euclidean norm from CAPD.	cout << "Phase 0 candidate ||x - P(x)|| = " << euclNorm((DDEs::Vector)segment - (DDEs::Vector)candidate) << endl;	capd::ddeshelper::plot_value("phase-0-candidate--", plot_dt, candidate, false);	capd::ddeshelper::plot_value("phase-0-Pcandidate--", plot_dt, segment, false);	setup.drawSolution("phase-0-full-trajectory--", trajectory);	setup.drawDelayMap("phase-0-plot-phasespace--", trajectory);	cout << "phase 1 start - refine by Newton method" << endl;	DDEs::Matrix DP(setup.M(), setup.M()), V(setup.M(), setup.M());	DDEs::Real diff = 1.0;	for (int z = 0; z < PARAMS::NUM_NEWTON; ++z){		diff = setup.refinePeriodic(cout, section, candidate, V, DP);		if (diff < PARAMS::EPS_NEWTON) break;	}	trajectory = setup.poincare(section, candidate, segment);  // segment = P(candidate)	cout << "Phase 1 candidate ||x - P(x)|| = " << euclNorm((DDEs::Vector)segment - (DDEs::Vector)candidate) << endl;	setup.drawSolution("phase-1-full-trajectory--", trajectory);	setup.drawDelayMap("phase-1-plot-phasespace--", trajectory);	cout << "Phase 2 start - finding good section and coordinates on it." << endl;	// Find a section - left eigenvector of V corresponding to eigenvalue 1	// TODO:	// Find section coordinates - recompute JacPoincare to the selected section,	// then and orthonormalize selected number of dominant-eigenvalue eigenvectors of DP	// TODO:	cout << "Phase 3 start - preparing initial sets for the rigorous part." << endl;	// TODO:	cout << "Done. Please run the following command to compute rigotrous inverse of the coordinate matrix:" << endl;	// TODO:	return candidate;}
