// see this file for the initialization of the library for this particular example
#include "setup.h"

using namespace std;
using namespace capd;

///////////////////////////////////////////
///////////////////////////////////////////
int main(){
	double delay = 1.0;
	double a = 3.5;
	double b = 7.95;
	int p = 32, n = 4;
	double h = delay / p;

	DS::Grid grid(h);
	DS::TimePoint tau = grid(p); 	// tau must be TimePoint from the grid! Not double!
	DS::TimePoint zero = grid(0);	// tau must be TimePoint from the grid! Not double!

	DS::Eq eq(a, b);
	DS::DDEq dde(eq, tau);
	DS::Solver dde_solver(dde);
	DS::Section S(1, 0, a/b);
	DS::PoincareMap P(dde_solver, S, capd::poincare::MinusPlus);
	P.setRequiredSteps(0);
	P.setMaxSteps(4*p);

	capd::DMap y_0("var:t;par:c,tau;fun:((t + tau)/tau)*c;", n+1);
	y_0.setParameter("c", 1.6*a/b);
	y_0.setParameter("tau", tau);
	DS::Solution X(-tau, zero, n, y_0);
	auto solution_A = X; dde_solver(solution_A, 4*tau);
	capd::ddeshelper::plot_value("periodic-candidate-A-", 0.1*h, solution_A);
	X = P(X);

	// We create a helper, that will allow to apply the Newton's algorithm easier
	// You can do this by yourself, but it is already implemented in the helper library.
	// Note: we don't require any steps (0), and we expect the return time to be at most 2
	//       the negative numbers means the full delays. This is just for convenience.
	DHelper helper({a, b, tau}, p, n, 0, -2);
	helper.setCrossingDirection(capd::poincare::MinusPlus);
	// this will hold the C^1 data - the Jacobian of the Poincare map at periodic point and the monodromy matrix
	capd::DMatrix V(X.dimension(), X.dimension()), DP = V;
	// this is just to hold extra debug output from refinePeriodic(), so it does not show on screen.
	std::ostringstream devnull;
	auto diff = helper.refinePeriodic(cerr, S, X, V, DP);
	// we refine with a much greater accuracy than in the first example
	double epsilon = 1e-10;
	while (std::abs(diff) > epsilon){
		cout << "diff (1) = " << diff << endl;
		diff = helper.refinePeriodic(cerr, S, X, V, DP);
	}
	cout << "Final refine diff = " << diff << endl;
	solution_A = X; dde_solver(solution_A, 4*tau);
	capd::ddeshelper::plot_value("periodic-candidate-A-", 0.1*h, solution_A);

	helper.setRequiredSteps(-(n+1));
	helper.setMaximumSteps(-2*n);
	diff = helper.refinePeriodic(cerr, S, X, V, DP);
	while (std::abs(diff) > epsilon){
		cout << "diff (2) = " << diff << endl;
		diff = helper.refinePeriodic(cerr, S, X, V, DP);
	}
	solution_A = X; dde_solver(solution_A, 4*tau);
	capd::ddeshelper::plot_value("periodic-candidate-B-", 0.1*h, solution_A);

	// finally, we convert the nonrigorous candidate into a rigorous representation
	// that we will use in rigorous computations. We save it in a binary format to
	// not loose the precision. NOTE: binary representation might not be portable between
	// different operational systems!
	capd::DVector X_as_DVector = X;    // we convert to a ,,raw'' representation, i.e. a collection of coefficients, not a function in some space
	capd::IVector IX { X_as_DVector }; // we convert to intervals
	// I usually output a short representation (of a few numbers) to compare if my operations are ok.
	// It is becouse the vectors are usualy very long and you cannot see anything in output.
	// The ddeshelper has a small function for shortening the output (up to a given number of characters).
	cout << capd::ddeshelper::slice(X_as_DVector, 100) << endl;
	cout << capd::ddeshelper::slice(IX, 100) << endl;
	// the filename is arbitrary, I tend to put .ivector for files in a human
	// readable format, and .ivector.bin for binary
	capd::ddeshelper::saveBinary("x0.ivector.bin", IX);

	// finally, we will use the matrices DP and V to compute approximate eigenvalues
	// and eigenvectors that will be helpful in proving periodic orbits.

	V.Transpose(); // we need V transposed to compute left eigenvectors (see FoCM papers for explanation)
	// the left eigenvector corresponding to eigenvalue 1. will be used as a section.
	auto M = IX.dimension();		// we use fact that this is exactly the dimensions of V, DP, and X as a vector.
	capd::DVector real(M), imag(M); 		// those will store eigenvalues
	capd::DMatrix vreal(M, M), vimag(M, M); // those will store eigenvectors as columns.
	capd::alglib::computeEigenvaluesAndEigenvectors(V, real, imag, vreal, vimag);
	int HOW_MUCH_PRINT = 30;
	cout << "EIGENVALUES OF THE MONODROMY MATRIX:" << endl;
	for (int i = 0; i < HOW_MUCH_PRINT; i++)
		cout << real[i] << " + " << imag[i] << "i " << (imag[i] == 0 ? capd::ddeshelper::slice(vreal.column(i)) : "") << endl;

	capd::DVector section = vreal.column(1);	// this should be the real eigenvector corresponding to eigenvalue 1
	section *= -1.0;
	// section.normalize();						// normalize section
	auto value_on_section = section * X_as_DVector; // the dot product of section and solution representation
	// this is the most general constructor of so called C^n_p-section (or (d,p,n)-section.
	// Basically, it is a section that is good for our representation of functions.
	// The value 1 is the d: dimension of the system (in this example 1 - scalar equation).
	DS::Section good_section(1, p, n, section, value_on_section);
	cout << "good_section(X): " << good_section(X) << " (should be close to 0.)" << endl;
	// we compute DP, but to the new section:
	// TODO: make helper better, so it is easier to use!
	auto cpyX = X;
	diff = helper.refinePeriodic(cerr, good_section, cpyX, V, DP);
	cout << "diff = " << diff << " should be close to 0, if not, switch the direction!" << endl;
	// now eigenvalues of DP. It should have eigenvalue ~0. instead of 1.
	capd::alglib::computeEigenvaluesAndEigenvectors(DP, real, imag, vreal, vimag);
	for (int i = 0; i < HOW_MUCH_PRINT; i++)
		cout << real[i] << " + " << imag[i] << "i" << (imag[i] == 0 ? capd::ddeshelper::slice(vreal.column(i)) : "") << endl;

	vreal.column(0); // this will correspond to the dominant eigenvalue
	vreal.column(1); // we might as well save this, if we need to subdivide in more dimensions.

//	IMatrix C(helper.M(), helper.M());
//	C.setToIdentity();
//	C.column(0) = section;
//	C.column(1) = vreal.column(0);
//	C.column(2) = vreal.column(1);
//
//	MpReal::setDefaultPrecision(512);
//	MpIMatrix mpC(helper.M(), helper.M());
//	to_mpi_matrix(C, mpC, 512);
//	MpIMatrix mpCinv = capd::matrixAlgorithms::inverseMatrix(mpC);


	return 0;

}
///////////////////////////////////////////
///////////////////////////////////////////

