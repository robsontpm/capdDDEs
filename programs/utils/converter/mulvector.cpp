#include <iostream>
#include <iomanip>
#include "common.h"

int main(int argc, char** argv){
	cout.precision(512);
	int precission = 512;
	int d = 1;
	std::string opt_help = "";
	if (argc >= 2){
		opt_help = argv[1];
	}
	if (argc < 6 || opt_help == "-h" || opt_help == "--h" || opt_help == "--help" || opt_help == "help" ){
		cout << "Program to multiply vectors, by scalars" << endl;
		cout << "usage: convvector [dimension] [interval|bin] [in_filepath] [r] [out_filepath]" << endl;
		return 0;
	}

	std::istringstream issdim(argv[1]); issdim >> d;
	std::string in_mode = argv[2];
	std::string in_filepath = argv[3];
	std::string r_string = argv[4];
	std::string out_filepath = argv[5];

	if ( in_mode != "interval" && in_mode != "bin" ){
		std::cout << "Bad in mode: " << in_mode << ", available options: interval|bin" << endl;
		return -1;
	}
	MIVector v(d);
	if (in_mode == "interval") {
		cout << "reading IVector START" << endl;
		ifstream f(in_filepath);
		f >> v;
		f.close();
		cout << "reading IVector DONE" << endl;
	}else if (in_mode == "bin") {
		cout << "reading BINARY IVector START" << endl;
		capd::ddeshelper::readBinary(in_filepath, v);
		cout << "reading BINARY IVector DONE" << endl;
	}

	interval r = 1.0;
	istringstream in_r(r_string);
	in_r >> r;

	auto out_v = r * v;

	if (in_mode == "interval") {
		cout << "Output IVector START" << endl;
		ofstream f(out_filepath); f.precision(16);
		f << out_v;
		f.close();
		cout << "Output IVector DONE" << endl;
	}else if (in_mode == "bin") {
		cout << "Output BINARY IVector START" << endl;
		capd::ddeshelper::saveBinary(out_filepath, out_v);
		cout << "Output BINARY IVector DONE" << endl;
	}
	cout << "# program " << argv[0] << " DONE" << endl;
	return 0;}
