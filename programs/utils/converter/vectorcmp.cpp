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
		cout << "Program to compare two ivectors. " << endl;
		cout << "usage: vectorcmp [dimension] [interval|bin] [in_prefix_A] [in_prefix_B]" << endl;
		return 0;
	}

	std::istringstream issdim(argv[1]); issdim >> d;
	std::string in_mode = argv[2];
	std::string in_path_A = argv[3];
	std::string in_path_B = argv[4];

	MIVector a(d), b(d);

	if ( in_mode != "interval" && in_mode != "bin"){
		std::cout << "Bad in mode: " << in_mode << ", available options: interval|bin" << endl;
		return -1;
	}

	if (in_mode == "interval") {
		cout << "reading plain text interval START" << endl;
		{ std::ifstream in(in_path_A); in >> a; in.close(); }
		{ std::ifstream in(in_path_B); in >> b; in.close(); }
		cout << "reading plain text interval DONE" << endl;
	}else if (in_mode == "bin") {
		cout << "reading BINARY START" << endl;
		capd::ddeshelper::readBinary(in_path_A, a);
		capd::ddeshelper::readBinary(in_path_B, b);
		cout << "reading BINARY DONE" << endl;
	}

	cout << "Comparison between:" << endl;
	cout << "    a: " << in_path_A << endl;
	cout << " and " << endl;
	cout << "    b: " << in_path_B << endl << endl;

	MIVector diffs = a - b; Interval sum = 0.0;
	for (int i = 0; i < d; ++i){
		Interval v = abs(diffs[i]); sum += v;
		cout << "at " << strfix(i, 5) << " |a - b| = " << printInterval(v) << endl;
	}
	cout << "Taxi norm a-b: " << printInterval(sum) << endl << endl;

	std::vector<int> counter; counter.resize(7, 0);
	for (int i = 0; i < d; ++i){
		auto relop = check_inclusion(a[i], b[i], counter);
		cout << "at " << i << " a rel b" << endl;
		cout << printInterval(a[i]) << endl;
		cout << relop.second << endl;
		cout << printInterval(b[i]) << endl;
		cout << endl;
	}
	print_counter(cout, counter);

	cout << "Taxi norm a-b: " << printInterval(sum) << endl << endl;

	return 0;}
