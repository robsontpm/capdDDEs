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
	if (argc < 5 || opt_help == "-h" || opt_help == "--h" || opt_help == "--help" || opt_help == "help" ){
		cout << "Program to check if two d x d matrices are good A A^{-1} pair. " << endl;
		cout << "usage: invmatrixtest [dimension] [interval|mp|bin|mpbin] [in_filepath_A] [in_filepath_A^{-1}] [precission(default: 512)]" << endl;
		return 0;
	}

	if (argc > 5){
		std::istringstream issprec(argv[5]); issprec >> precission;
	}

	std::istringstream issdim(argv[1]); issdim >> d;
	std::string in_mode = argv[2];
	std::string in_filepath_A = argv[3];
	std::string in_filepath_Ainv = argv[4];

	if ( in_mode != "interval" && in_mode != "mp" && in_mode != "bin"  && in_mode != "mpbin"){
		std::cout << "Bad in mode: " << in_mode << ", available options: interval|mp|bin|mpbin" << endl;
		return -1;
	}

	if (in_mode == "interval") {
		cout << "reading IMatrix START" << endl;
		MIMatrix iA(d, d);
		MIMatrix iAinv(d, d);
		{ ifstream f(in_filepath_A); f >> iA; f.close(); }
		{ ifstream f(in_filepath_Ainv); f >> iAinv; f.close(); }
		testInverse(iA, iAinv);
		cout << "reading IMatrix DONE" << endl;
	}else if (in_mode == "mp") {
		cout << "reading MpIMatrix START" << endl;
		int old_precision = MpReal::getDefaultPrecision();
		MpReal::setDefaultPrecision(precission);
		MpIMatrix iA(d, d);
		MpIMatrix iAinv(d, d);
		{ ifstream f(in_filepath_A); f >> iA; f.close(); }
		{ ifstream f(in_filepath_Ainv); f >> iAinv; f.close(); }
		testInverse(iA, iAinv);
		MpReal::setDefaultPrecision(old_precision);
		cout << "reading MpIMatrix DONE" << endl;
	}else if (in_mode == "bin") {
		cout << "reading BINARY IMatrix START" << endl;
		MIMatrix iA(d, d);
		MIMatrix iAinv(d, d);
		capd::ddeshelper::readBinary(in_filepath_A, iA);
		capd::ddeshelper::readBinary(in_filepath_Ainv, iAinv);
		testInverse(iA, iAinv);
		cout << "reading BINARY IMatrix DONE" << endl;
	}else if (in_mode == "mpbin") {
		cout << "reading BINARY MpIMatrix START" << endl;
		int old_precision = MpReal::getDefaultPrecision();
		MpReal::setDefaultPrecision(precission);
		MpIMatrix iA(d, d);
		MpIMatrix iAinv(d, d);
		capd::ddeshelper::readBinary(in_filepath_A, iA);
		capd::ddeshelper::readBinary(in_filepath_Ainv, iAinv);
		MpReal::setDefaultPrecision(old_precision);
		cout << "reading BINARY MpIMatrix DONE" << endl;
	}
	cout << "# program " << argv[0] << " DONE" << endl;
	return 0;}
