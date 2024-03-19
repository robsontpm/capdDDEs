#ifndef EXAMPLES_ROSSLER_ODE_VS_DDE_CODE_GNUPLOT_H_
#define EXAMPLES_ROSSLER_ODE_VS_DDE_CODE_GNUPLOT_H_

#include "setup.h"

std::string gnuplotHeader(double box[], int slices, std::string const& out_filename);

// This is just to plot some pictures
template<typename VectorType>
std::string gnuplotHullBox(
			VectorType const& A,
			int slice, int n_slices,
			std::string const& comment){
	std::ostringstream oss; oss.precision(15);

	const int nominally_unstable_index = 1;

	double x1 = A[1].leftBound();
	double x2 = A[1].rightBound();

	double y1 = A[2].leftBound();
	double y2 = A[2].rightBound();

	oss << "\n";
	oss << "# " << comment << " slice " << (slice >= n_slices ? slice - n_slices : slice) + 1 << " of " << n_slices << "\n";
	oss << "set obj rect from " << x1 << "," << y1 << " to " << x2 << "," << y2 << " fs solid 1.0 fc palette cb " << slice << "\n";


	return oss.str();
}

template<typename Vector1, typename Vector2>
void plotTrappingRegion(double box[], std::string name, std::vector<Vector1> const& x, std::vector<Vector2> const& Px){
	std::ofstream gp(name + ".gp");
	gp.precision(16);
	gp << gnuplotHeader(box, x.size(), name + ".png");
	for (int i = 0; i < x.size(); ++i)
		gp << gnuplotHullBox(x[i], i + x.size(), x.size(), "") << "\n";

	for (int i = 0; i < Px.size(); ++i)
		gp << gnuplotHullBox(Px[i], i, Px.size(), "");

	gp << "plot 0" << std::endl;
	gp.close();

	#ifdef DDES_ALLOW_SYSTEM
	capd::ddeshelper::runSystemCommand("gnuplot " + name + ".gp");
	#endif
}

#endif // EXAMPLES_ROSSLER_ODE_VS_DDE_CODE_GNUPLOT_H_
