#define WITH_GNUPLOT

#include <iostream>
#include <iomanip>
#include <capd/capdlib.h>
#include <capd/ddes/ddeslib.h>
#include <capd/ddeshelper/ddeshelperlib.h>
#include <capd/ddeshelper/DDEHelperRigorous.hpp>
#include "constants.h"
#include <algorithm>

// TODO: make it work later

using namespace std;

typedef capd::intervals::Interval<double, capd::rounding::DoubleRounding>  Interval;
typedef capd::ddes::MackeyGlass<Interval, Interval> Eq;
typedef capd::ddeshelper::RigorousHelper<Eq> Setup;

///////////////////////////////////////////
///////////////////////////////////////////
int main(int argc, char** argv){

	Setup::Vector test_r0({0.0, 0.1, 0.1, 1.0, 1.0});
	capd::ddeshelper::CubeSet<Setup::Vector> cubes(test_r0);
	cubes.insert({0,0,0});
	cubes.insert({0,1,0});
	cubes.insert({0,-1,0});
	cubes.insert({0,0,-1});
	cubes.insert({0,0,1});
	for (auto it = cubes.csbegin(); it != cubes.csend(); ++it){
		std::cout << (*it) << std::endl;
	}
	std::cout << "All cube signatures shown!" << endl;
	for (auto it = cubes.begin(); it != cubes.end(); ++it){
		std::cout << (*it) << " " << Setup::Vector(it) << std::endl;
	}
	std::cout << "All cubes shown!" << endl;
	for (auto& it: cubes){
		std::cout << it << " " << cubes.get_dx(it) << std::endl;
	}
	std::cout << "All cubes shown!" << endl;
	std::cout << "Insert BOX test (good box)" << endl;
	Setup::Vector v(5);
	capd::ddeshelper::CubeSet<Setup::Vector> cubes2(test_r0);
	v[0] = Setup::Scalar(-0.2, 0.2);
	v[1] = Setup::Scalar(-0.3, 0.3);
	v[2] = Setup::Scalar(-0.2, 0.2);
	v[3] = Setup::Scalar(-0.03, 0.03);
	v[4] = Setup::Scalar(-0.5, 0.5);
	auto last_cut_needed = cubes2.insert_cover(v, 4);
	cout << "last_cut_needed = " << last_cut_needed << ", should be smaller than 4" << std::endl;
	for (auto it = cubes2.csbegin(); it != cubes2.csend(); ++it){
		std::cout << (*it) << std::endl;
	}
	std::cout << "All cube signatures shown!" << endl;
	for (auto it = cubes2.begin(); it != cubes2.end(); ++it){
		std::cout << (*it) << " " << Setup::Vector(it) << std::endl;
	}
	std::cout << "All cubes shown!" << endl;
	for (auto& it: cubes2){
		std::cout << it << " " << cubes2.get_dx(it) << std::endl;
	}
	std::cout << "All cubes shown!" << endl;
	std::cout << "Insert BOX test (too big box)" << endl;
	capd::ddeshelper::CubeSet<Setup::Vector> cubes3(test_r0);
	v[3] = Setup::Scalar(-0.1, 5.1);
	v[4] = Setup::Scalar(-1.5, 0.5);
	last_cut_needed = cubes3.insert_cover(v, 4);
	cout << "last_cut_needed = " << last_cut_needed << ", will be >= 4" << std::endl;
	for (auto it = cubes3.csbegin(); it != cubes3.csend(); ++it){
		std::cout << (*it) << std::endl;
	}
	std::cout << "All cube signatures shown!" << endl;
	for (auto it = cubes3.begin(); it != cubes3.end(); ++it){
		std::cout << (*it) << " " << Setup::Vector(it) << std::endl;
	}
	std::cout << "All cubes shown!" << endl;
	for (auto& it: cubes3){
		std::cout << it << " " << cubes3.get_dx(it) << std::endl;
	}

	std::cout << "The same cube sets? (cubes3 == cubes3)" << (cubes3 == cubes3) << std::endl;
	std::cout << "The same cube sets? (cubes2 == cubes2)" << (cubes2 == cubes2) << std::endl;
	std::cout << "The same cube sets? (cubes3 == cubes2)" << (cubes3 == cubes2) << std::endl;

	std::cout << "Out operator test" << endl;
	std::cout << "cubes: " << cubes << endl;
	capd::ddeshelper::CubeSet<Setup::Vector> cubes_empty(test_r0);
	std::cout << "cubes_empty: " << cubes_empty << endl;
	std::cout << "Out operator test DONE" << endl;

	std::cout << "In operator test" << endl;
	std::istringstream no_ws_input("{[0,-1,0],[0,0,-1],[0,0,0],[0,0,1],[0,1,0]} asadasd wsdaw wada");
	std::istringstream ws_input("{    \n   [0,-1,0]\n,[0,0,-1]\n     ,[0,0,0]     \n,[0,0,1]   \n    ,\t[0,1,0]  \t\n\n}\n\n wdadwa  w ada");
	capd::ddeshelper::CubeSet<Setup::Vector> cubes4(test_r0);
	no_ws_input >> cubes4;
	capd::ddeshelper::CubeSet<Setup::Vector> cubes5(test_r0);
	ws_input >> cubes5;
	cout << "cubes4: " << cubes4 << endl;
	cout << "cubes5: " << cubes5 << endl;
	std::cout << "The same cube sets? (cubes4 == cubes5)" << (cubes4 == cubes5) << std::endl;
	std::cout << "The same cube sets? (cubes == cubes4)" << (cubes == cubes4) << std::endl;
	std::cout << "The same cube sets? (cubes == cubes5)" << (cubes == cubes5) << std::endl;
	std::ostringstream osstest;
	osstest << cubes5 << "\n";
	std::istringstream isstest(osstest.str());
	capd::ddeshelper::CubeSet<Setup::Vector> cubes6(test_r0);
	isstest >> cubes6;
	std::cout << "The same cube sets? (cubes6 == cubes5)" << (cubes6 == cubes5) << std::endl;
	std::cout << "In operator test DONE" << endl;

	return 0;

	return 0;
}
