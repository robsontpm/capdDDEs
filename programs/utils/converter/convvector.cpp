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
	if (argc < 7 || opt_help == "-h" || opt_help == "--h" || opt_help == "--help" || opt_help == "help" ){
		cout << "Program to convert vectors to various data types. " << endl;
		cout << "usage: convvector [dimension] [double|interval|mp|bin|mpbin] [in_filepath] [conv|makeradius] [interval|mp|bin|mpbin] [out_filepath] [precission(default: 512)]" << endl;
		return 0;
	}

	if (argc > 7){
		std::istringstream issprec(argv[7]); issprec >> precission;
	}

	std::istringstream issdim(argv[1]); issdim >> d;
	std::string in_mode = argv[2];
	std::string in_filepath = argv[3];
	std::string op_mode = argv[4];
	std::string out_mode = argv[5];
	std::string out_filepath = argv[6];

	if ( in_mode != "double" && in_mode != "interval" && in_mode != "mp" && in_mode != "bin"  && in_mode != "mpbin"){
		std::cout << "Bad in mode: " << in_mode << ", available options: double|interval|mp|bin|mpbin" << endl;
		return -1;
	}
	if ( out_mode != "interval" && out_mode != "mp" && out_mode != "bin"  && out_mode != "mpbin"){
		std::cout << "Bad out mode: " << in_mode << ", available options: interval|mp|bin|mpbin" << endl;
		return -1;
	}
	if ( op_mode != "conv" && op_mode != "makeradius"){
		std::cout << "Bad operation: " << op_mode << ", available options: conv|makeradius (maybe add future operations!)" << endl;
		return -1;
	}

	MpIVector v(d);
	if (in_mode == "double"){
		cout << "reading DVector START" << endl;
		MDVector dv(d);
		ifstream f(in_filepath);
		f >> dv;
		to_mpi_vector(dv, v, precission);
		cout << "reading DVector DONE" << endl;
	}else if (in_mode == "interval") {
		cout << "reading IVector START" << endl;
		MIVector iv(d);
		ifstream f(in_filepath);
		f >> iv;
		to_mpi_vector(iv, v, precission);
		f.close();
		cout << "reading IVector DONE" << endl;
	}else if (in_mode == "mp") {
		cout << "reading MpIVector START" << endl;
		ifstream f(in_filepath);
		f >> v;
		f.close();
		cout << "reading MpIVector DONE" << endl;
	}else if (in_mode == "bin") {
		cout << "reading BINARY IVector START" << endl;
		MIVector iv(d);
		capd::ddeshelper::readBinary(in_filepath, iv);
		to_mpi_vector(iv, v, precission);
		cout << "reading BINARY IVector DONE" << endl;
	}else if (in_mode == "mpbin") {
		cout << "reading BINARY MpIVector START" << endl;
		capd::ddeshelper::readBinary(in_filepath, v);
		cout << "reading BINARY MpIVector DONE" << endl;
	}

	MpIVector out_v(d);
	if (op_mode == "makeradius"){
		cout << "Operation: makeradius START" << endl;
		out_v = MpInterval(-1.0, 1.0) * v;
		cout << "Operation: makeradius DONE" << endl;
	} else {
		cout << "Operation: just convert START" << endl;
		out_v = v;
		cout << "Operation: just convert DONE" << endl;
	}

	if (out_mode == "interval") {
		cout << "Output IVector START" << endl;
		MIVector iv(d);
		from_mpi_vector(out_v, iv);
		ofstream f(out_filepath); f.precision(16);
		f << iv;
		f.close();
		cout << "Output IVector DONE" << endl;
	}else if (out_mode == "mp") {
		cout << "Output MpIVector START" << endl;
		ofstream f(out_filepath); f.precision(precission);
		f << out_v;
		f.close();
		cout << "Output MpIVector DONE" << endl;
	}else if (out_mode == "bin") {
		cout << "Output BINARY IVector START" << endl;
		MIVector iv(d);
		from_mpi_vector(out_v, iv);
		capd::ddeshelper::saveBinary(out_filepath, iv);
		cout << "Output BINARY IVector DONE" << endl;
	}else if (out_mode == "mpbin") {
		cout << "Output BINARY MpIVector START" << endl;
		capd::ddeshelper::saveBinary(out_filepath, out_v);
		cout << "Output BINARY MpIVector DONE" << endl;
	}
	cout << "# program " << argv[0] << " DONE" << endl;
	return 0;}
