#include <iostream>
#include <iomanip>
#include "common.h"

// TODO: TU SKONCZYLEM

int main(int argc, char** argv){
	cout.precision(15);
	int d = 1;
	std::string opt_help = "";
	if (argc >= 2){
		opt_help = argv[1];
	}
	if (argc < 6 || opt_help == "-h" || opt_help == "--h" || opt_help == "--help" || opt_help == "help" ){
		cout << "Program to convert r vector to a matrix R and R^-1, with r[i] and r[i]^{-1} on diagonal" << endl;
		cout << "r is symetrized and max is taken before making R. " << endl;
		cout << "usage: rmatrix [dimension] [interval|bin] [in_r_path] [out_R_path] [in_invR_path] [start_at]" << endl;
		return 0;
	}

	std::istringstream issdim(argv[1]); issdim >> d;
	std::string mode = argv[2];

	std::string in_r_path = argv[3];

	std::string out_R_path = argv[4];
	std::string out_invR_path = argv[5];

	int start_at = 0;
	if (argc > 6) { std::istringstream iss(argv[6]); iss >> start_at;}

	if ( mode != "interval" && mode != "bin"){
		std::cout << "Bad in mode: " << mode << ", available options: interval|bin" << endl;
		return -1;
	}

	MIMatrix R(d, d), invR(d, d);
	MIVector in_r(d);
	if (mode == "interval") {
		cout << "reading plain text interval START" << endl;
		{ std::ifstream in(in_r_path); in >> in_r; in.close(); }
		cout << "reading plain text interval DONE" << endl;
	}else if (mode == "bin") {
		cout << "reading BINARY START" << endl;
		capd::ddeshelper::readBinary(in_r_path, in_r);
		cout << "reading BINARY DONE" << endl;
	}

	in_r = capd::vectalg::abs(in_r);
	R.setToIdentity(); invR.setToIdentity();
	for (int i = start_at; i < d; ++i){
		if (in_r[i].rightBound() > 0.){
			R[i][i] = in_r[i].rightBound();
			invR[i][i] = MIVector::ScalarType(1.0) / MIVector::ScalarType(R[i][i]);
		}
	}

	if (mode == "interval") {
		cout << "writing plain text interval START" << endl;
		{ std::ofstream out(out_R_path); 	out.precision(16); out << R; 	out.close(); }
		{ std::ofstream out(out_invR_path); out.precision(16); out << invR; out.close(); }
		cout << "writing plain text interval DONE" << endl;
	}else if (mode == "bin") {
		cout << "reading BINARY START" << endl;
		capd::ddeshelper::saveBinary(out_R_path, R);
		capd::ddeshelper::saveBinary(out_invR_path, invR);
		cout << "writing BINARY DONE" << endl;
	}

	cout << "# program " << argv[0] << " DONE" << endl;

	return 0;}
