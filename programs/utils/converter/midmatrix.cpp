#include <iostream>
#include <iomanip>
#include "common.h"

int main(int argc, char** argv){
	cout.precision(16);
	int d = 1;
	std::string opt_help = "";
	if (argc >= 2){
		opt_help = argv[1];
	}
	if (argc < 5 || opt_help == "-h" || opt_help == "--h" || opt_help == "--help" || opt_help == "help" ){
		cout << "Program to extract size d x d middle matrix and save as double. " << endl;
		cout << "usage: midmatrix [dimension] [interval|bin] [in_filepath] [out_filepath]" << endl;
		return 0;
	}

	std::istringstream issdim(argv[1]); issdim >> d;
	std::string in_mode = argv[2];
	std::string in_filepath = argv[3];
	std::string out_filepath = argv[4];

	if (in_mode != "interval" && in_mode != "bin"){
		std::cout << "Bad in mode: " << in_mode << ", available options: interval|bin" << endl;
		return -1;
	}

	MIMatrix A(d, d);
	if (in_mode == "interval") {
		cout << "reading IMatrix START" << endl;
		ifstream f(in_filepath);
		f >> A;
		f.close();
		cout << "reading IMatrix DONE" << endl;
	}else if (in_mode == "bin") {
		cout << "reading BINARY IMatrix START" << endl;
		capd::ddeshelper::readBinary(in_filepath, A);
		cout << "reading BINARY IMatrix DONE" << endl;
	}

	MDMatrix outA(d, d);
	for (int i = 0; i < d; ++i)
		for (int j = 0; j < d; ++j)
			outA[i][j] = (A[i][j].leftBound() + A[i][j].rightBound()) / 2.0;

	cout << "Output DMatrix START" << endl;
	ofstream f(out_filepath); f.precision(16);
	f << outA;
	f.close();
	cout << "Output DMatrix DONE" << endl;
	cout << "# program " << argv[0] << " DONE" << endl;

	return 0;
}
