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
		cout << "Program to convert and eventually manipulate matrices of size d x d (inversion with multiprecission!). " << endl;
		cout << "usage: convmatrix [dimension] [double|interval|mp|bin|mpbin] [in_filepath] [conv|inv|T] [double|interval|mp|bin|mpbin] [out_filepath] [precission(default: 512)]" << endl;
		cout << "WARNING: conversion to double available only from interval and bin!" << endl;
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
	if ( out_mode != "double" && out_mode != "interval" && out_mode != "mp" && out_mode != "bin"  && out_mode != "mpbin"){
		std::cout << "Bad out mode: " << in_mode << ", available options: double|interval|mp|bin|mpbin" << endl;
		return -1;
	}
	if ( op_mode != "inv" && op_mode != "conv" && op_mode != "T"  && op_mode != "transpose"){
		std::cout << "Bad operation: " << op_mode << ", available options: conv|inv|(T|transpose)" << endl;
		return -1;
	}

	if ( out_mode == "double" && !(op_mode == "conv" && (in_mode == "bin" || in_mode == "interval" ))){
		std::cout << "Only conversion operation available to double and only from: bin or interval!" << endl;
		return -1;
	}

	MIMatrix IA(d, d); // this is a dirty fast hack to not alter he other code. In fact, I should be able to pull iA outside, TODO: refactor
	MpIMatrix A(d, d);
	if (in_mode == "double"){
		cout << "reading DMatrix START" << endl;
		MDMatrix dA(d, d);
		ifstream f(in_filepath);
		f >> dA;
		to_mpi_matrix(dA, A, precission);
		cout << "reading DMatrix DONE" << endl;
	}else if (in_mode == "interval") {
		cout << "reading IMatrix START" << endl;
		MIMatrix iA(d, d);
		ifstream f(in_filepath);
		f >> iA; IA = iA;
		to_mpi_matrix(iA, A, precission);
		f.close();
		cout << "reading IMatrix DONE" << endl;
	}else if (in_mode == "mp") {
		cout << "reading MpIMatrix START" << endl;
		ifstream f(in_filepath);
		f >> A;
		f.close();
		cout << "reading MpIMatrix DONE" << endl;
	}else if (in_mode == "bin") {
		cout << "reading BINARY IMatrix START" << endl;
		MIMatrix iA(d, d);
		capd::ddeshelper::readBinary(in_filepath, iA);
		IA = iA;
		to_mpi_matrix(iA, A, precission);
		cout << "reading BINARY IMatrix DONE" << endl;
	}else if (in_mode == "mpbin") {
		cout << "reading BINARY MpIMatrix START" << endl;
		capd::ddeshelper::readBinary(in_filepath, A);
		cout << "reading BINARY MpIMatrix DONE" << endl;
	}

	MpIMatrix outA(d, d);
	if (op_mode == "T" || op_mode == "transpose"){
		cout << "Operation: transpose START" << endl;
		outA = A;
		outA.transpose();
		cout << "Operation: transpose DONE" << endl;
	}else if (op_mode == "inv"){
		cout << "Operation: inverse START" << endl;
		outA = mpi_inverse(A);
		testInverse(A, outA);
		cout << "Operation: inverse DONE" << endl;
	}else{
		cout << "Operation: just convert START" << endl;
		outA = A;
		cout << "Operation: just convert DONE" << endl;
	}

	if (out_mode == "interval") {
		cout << "Output IMatrix START" << endl;
		MIMatrix iA(d, d);
		from_mpi_matrix(outA, iA);
		ofstream f(out_filepath); f.precision(16);
		f << iA;
		f.close();
		cout << "Output IMatrix DONE" << endl;
	}else if (out_mode == "mp") {
		cout << "Output MpIMatrix START" << endl;
		ofstream f(out_filepath); f.precision(precission);
		f << outA;
		f.close();
		cout << "Output MpIMatrix DONE" << endl;
	}else if (out_mode == "bin") {
		cout << "Output BINARY IMatrix START" << endl;
		MIMatrix iA(d, d);
		from_mpi_matrix(outA, iA);
		capd::ddeshelper::saveBinary(out_filepath, iA);
		cout << "Output BINARY IMatrix DONE" << endl;
	}else if (out_mode == "mpbin") {
		cout << "Output BINARY MpIMatrix START" << endl;
		capd::ddeshelper::saveBinary(out_filepath, outA);
		cout << "Output BINARY MpIMatrix DONE" << endl;
	}else if (out_mode == "double") {
		cout << "Output DMatrix START" << endl;
		auto midIA = midMatrix(IA);
		MDMatrix dA(d, d);
		auto ia = midIA.begin();
		auto id = dA.begin();
		for (; ia != midIA.end(); ++ia, ++id)
			(*id) = (*ia).mid().leftBound();
		ofstream f(out_filepath); f.precision(16);
		f << dA;
		f.close();
		cout << "Output DMatrix DONE" << endl;
	}

	cout << "# program " << argv[0] << " DONE" << endl;
	return 0;}
