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
	if (argc < 6 || opt_help == "-h" || opt_help == "--h" || opt_help == "--help" || opt_help == "help" ){
		cout << "Program to make a hull of ivectors. " << endl;
		cout << "usage: vectorhull [dimension] [interval|bin] [out_prefix] [in_prefix_1] [...]" << endl;
		return 0;
	}

	std::istringstream issdim(argv[1]); issdim >> d;
	std::string mode = argv[2];
	std::vector<std::string> in_paths;
	std::string out_path = argv[3];

	if ( mode != "interval" && mode != "bin"){
		std::cout << "Bad in mode: " << mode << ", available options: interval|bin" << endl;
		return -1;
	}

	for (int i = 4; i < argc; ++i)
		in_paths.push_back(std::string(argv[i]));

	MIVector zero(d);
	std::vector<MIVector> in_vectors(in_paths.size(), zero);
	for (int i = 0; i < in_paths.size(); i++){
		if (mode == "interval") {
			cout << "reading plain text interval " << i << " START" << endl;
			{ std::ifstream in(in_paths[i]); in >> in_vectors[i]; in.close(); }
			cout << "reading plain text interval " << i << " DONE" << endl;
		}else if (mode == "bin") {
			cout << "reading BINARY " << i << " START" << endl;
			capd::ddeshelper::readBinary(in_paths[i], in_vectors[i]);
			cout << "reading BINARY " << i << " DONE" << endl;
		}
	}

	MIVector hull = in_vectors[0];
	for (int i = 1; i < in_vectors.size(); ++i)
		hull = capd::vectalg::intervalHull(hull, in_vectors[i]);

	if (mode == "interval") {
		cout << "writing plain text interval START" << endl;
		{ std::ofstream out(out_path); out.precision(16); out << hull; out.close(); }
		cout << "writing plain text interval DONE" << endl;
	}else if (mode == "bin") {
		cout << "reading BINARY START" << endl;
		capd::ddeshelper::saveBinary(out_path, hull);
		cout << "writing BINARY DONE" << endl;
	}
	cout << "# program " << argv[0] << " DONE" << endl;
	return 0;}
