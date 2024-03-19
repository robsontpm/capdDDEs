#include "gnuplot.h"

std::string gnuplotHeader(double box[], int slices, std::string const& out_filename){
	std::ostringstream oss;

	oss << "set terminal png size 2400,1600" << "\n";
	oss << "set output '" << out_filename << "'" << "\n";
	oss << "set xrange [" << (box[0] - 1.) << ":" << (box[1] + 1.) << "]" << "\n";
	oss << "set yrange [" << box[2] - abs(box[2]) << ":" << box[3] + abs(box[3]) << "]" << "\n";
	oss << "set palette model RGB defined ( 0 'yellow', 1 'red', 1 'blue', 2 'green')  " << "\n";
	oss << "set cbrange [0:" << 2*slices << "]" << "\n";

	return oss.str();
}
