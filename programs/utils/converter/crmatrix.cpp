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
	if (argc < 8 || opt_help == "-h" || opt_help == "--h" || opt_help == "--help" || opt_help == "help" ){
		cout << "Program to convert C * r0 (with C^{-1}), to (CR) * [-1,1]^d, with (R^{-1} C^{-1})." << endl;
		cout << "usage: crmatrix [dimension] [interval|bin] [in_C_path] [in_invC_path] [in_r0_path] [out_CR_path] [in_invRinvC_path] [start_at]" << endl;
		return 0;
	}

	std::istringstream issdim(argv[1]); issdim >> d;
	std::string mode = argv[2];

	std::string in_C_path = argv[3];
	std::string in_invC_path = argv[4];
	std::string in_r0_path = argv[5];

	std::string out_CR_path = argv[6];
	std::string out_invRinvC_path = argv[7];

	int start_at = 0;
	if (argc > 8) { std::istringstream iss(argv[8]); iss >> start_at; }

	if ( mode != "interval" && mode != "bin"){
		std::cout << "Bad in mode: " << mode << ", available options: interval|bin" << endl;
		return -1;
	}

	MIMatrix in_C(d, d), in_invC(d, d);
	MIMatrix CR(d, d), invRinvC(d, d);
	MIVector in_r0(d);
	MIMatrix R(d, d), invR(d, d);
	if (mode == "interval") {
		cout << "reading plain text interval START" << endl;
		{ std::ifstream in(in_C_path); in >> in_C; in.close(); }
		{ std::ifstream in(in_invC_path); in >> in_invC; in.close(); }
		{ std::ifstream in(in_r0_path); in >> in_r0; in.close(); }
		cout << "reading plain text interval DONE" << endl;
	}else if (mode == "bin") {
		cout << "reading BINARY START" << endl;
		capd::ddeshelper::readBinary(in_C_path, in_C);
		capd::ddeshelper::readBinary(in_invC_path, in_invC);
		capd::ddeshelper::readBinary(in_r0_path, in_r0);
		cout << "reading BINARY DONE" << endl;
	}

	in_r0 = capd::vectalg::abs(in_r0);
	R.setToIdentity(); invR.setToIdentity();
	for (int i = start_at; i < d; ++i){
		if (in_r0[i].rightBound() > 0.){
			R[i][i] = in_r0[i].rightBound();
			invR[i][i] = MIVector::ScalarType(1.0) / MIVector::ScalarType(R[i][i]);
		}
	}

	CR = in_C * R;
	invRinvC = invR * in_invC;

	if (mode == "interval") {
		cout << "writing plain text interval START" << endl;
		{ std::ofstream out(out_CR_path); out.precision(16); out << CR; out.close(); }
		{ std::ofstream out(out_invRinvC_path); out.precision(16); out << invRinvC; out.close(); }
		cout << "writing plain text interval DONE" << endl;
	}else if (mode == "bin") {
		cout << "reading BINARY START" << endl;
		capd::ddeshelper::saveBinary(out_CR_path, CR);
		capd::ddeshelper::saveBinary(out_invRinvC_path, invRinvC);
		cout << "writing BINARY DONE" << endl;
	}
	cout << "# program " << argv[0] << " DONE" << endl;
	return 0;}
