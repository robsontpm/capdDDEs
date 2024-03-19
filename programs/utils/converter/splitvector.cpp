#include <iostream>
#include <iomanip>
#include "common.h"

int main(int argc, char** argv){
	int d = 1;
	std::string opt_help = "";
	if (argc >= 2){
		opt_help = argv[1];
	}
	if (argc < 6 || opt_help == "-h" || opt_help == "--h" || opt_help == "--help" || opt_help == "help" ){
		cout << "Program to split vector into ref + r, where ref is point vector and r is 0-centered. " << endl;
		cout << "usage: splitvector [dimension] [interval|bin] [in_filepath] [out_ref_filepath] [out_r_filepath]" << endl;
		return 0;
	}

	std::istringstream issdim(argv[1]); issdim >> d;
	std::string in_mode = argv[2];
	std::string in_filepath = argv[3];
	std::string out_ref_filepath = argv[4];
	std::string out_r_filepath = argv[5];

	if ( in_mode != "interval" && in_mode != "bin"){
		std::cout << "Bad in mode: " << in_mode << ", available options: interval|bin" << endl;
		return -1;
	}

	MIVector iv(d), r(d);
	if (in_mode == "interval") {
		cout << "reading IVector START" << endl;
		ifstream f(in_filepath);
		f >> iv;
		f.close();
		cout << "reading IVector DONE" << endl;
	}else if (in_mode == "bin") {
		cout << "reading BINARY IVector START" << endl;
		capd::ddeshelper::readBinary(in_filepath, iv);
		cout << "reading BINARY IVector DONE" << endl;
	}

	capd::vectalg::split(iv, r);

	if (in_mode == "interval") {
		cout << "Output IVector START" << endl;
		{ ofstream f(out_ref_filepath); f.precision(16); f << iv; f.close(); }
		{ ofstream f(out_r_filepath);   f.precision(16); f << r;  f.close(); }
		cout << "Output IVector DONE" << endl;
	}else if (in_mode == "bin") {
		cout << "Output BINARY IVector START" << endl;
		capd::ddeshelper::saveBinary(out_ref_filepath, iv);
		capd::ddeshelper::saveBinary(out_r_filepath, r);
		cout << "Output BINARY IVector DONE" << endl;
	}

	cout << "# program " << argv[0] << " DONE" << endl;}
