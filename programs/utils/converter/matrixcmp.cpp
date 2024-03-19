#include <iostream>
#include <iomanip>
#include "common.h"

int main(int argc, char** argv){
	cout.precision(15);
	int d = 1;
	std::string opt_help = "";
	if (argc >= 2){
		opt_help = argv[1];
	}
	if (argc < 5 || opt_help == "-h" || opt_help == "--h" || opt_help == "--help" || opt_help == "help" ){
		cout << "Program to compare two imatrices. " << endl;
		cout << "usage: vectorcmp [dimension] [interval|bin] [in_path_A] [in_path_B]" << endl;
		return 0;
	}

	std::istringstream issdim(argv[1]); issdim >> d;
	std::string in_mode = argv[2];
	std::string in_path_A = argv[3];
	std::string in_path_B = argv[4];

	MIMatrix A(d, d), B(d, d);

	if ( in_mode != "interval" && in_mode != "bin"){
		std::cout << "Bad in mode: " << in_mode << ", available options: interval|bin" << endl;
		return -1;
	}

	if (in_mode == "interval") {
		cout << "reading plain text interval START" << endl;
		{ std::ifstream in(in_path_A); in >> A; in.close(); }
		{ std::ifstream in(in_path_B); in >> B; in.close(); }
		cout << "reading plain text interval DONE" << endl;
	}else if (in_mode == "bin") {
		cout << "reading BINARY START" << endl;
		capd::ddeshelper::readBinary(in_path_A, A);
		capd::ddeshelper::readBinary(in_path_B, B);
		cout << "reading BINARY DONE" << endl;
	}

	cout << "Comparison between:" << endl;
	cout << "    A: " << in_path_A << endl;
	cout << " and " << endl;
	cout << "    B: " << in_path_B << endl << endl;

	std::vector<int> counter; counter.resize(7, 0);
	std::vector<Interval> col_sums;
	for (int col = 0; col < d; ++col){
		cout << "column " << col << "start" << endl;

		MIVector diffs = A.column(col) - B.column(col); Interval col_sum = 0.0;
		for (int i = 0; i < d; ++i){
			Interval v = abs(diffs[i]); col_sum += v;
			cout << "at " << strfix(i, 5) << " |a - b| = " << printInterval(v) << endl;
		}
		cout << "Taxi norm A.col(" << strfix(col, 4) << ") - B.col(" << strfix(col, 4) << ") << : " << printInterval(col_sum) << endl << endl;
		col_sums.push_back(col_sum);

		for (int i = 0; i < d; ++i){
			auto relop = check_inclusion(A[i][col], B[i][col], counter);
			cout << "column " << col << " at " << i << " A rel B" << endl;
			cout << printInterval(A[i][col]) << endl;
			cout << relop.second << endl;
			cout << printInterval(B[i][col]) << endl;
			cout << endl;
		}
		cout << "Taxi norm A.col(" << strfix(col, 4) << ") - B.col(" << strfix(col, 4) << ") << : " << printInterval(col_sum) << endl << endl;
		cout << "column " << col << "end" << endl << endl;
	}

	for (int col = 0; col < d; col++){
		cout << "Taxi norm A.col(" << strfix(col, 4) << ") - B.col(" << strfix(col, 4) << ") << : " << printInterval(col_sums[col]) << endl;
	}
	print_counter(cout, counter);
	cout << "# program " << argv[0] << " DONE" << endl;
	return 0;}
