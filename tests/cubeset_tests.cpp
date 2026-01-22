#define BOOST_TEST_MODULE CubeSetTests
#include <boost/test/included/unit_test.hpp>

#include <iostream>
#include <iomanip>
#include <capd/capdlib.h>
#include <capd/ddes/ddeslib.h>
#include <capd/ddeshelper/ddeshelperlib.h>
#include <capd/ddeshelper/DDEHelperRigorous.hpp>
#include <algorithm>

using namespace std;

typedef capd::intervals::Interval<double, capd::rounding::DoubleRounding>  Interval;
typedef capd::ddes::MackeyGlass<Interval, Interval> Eq;
typedef capd::ddeshelper::RigorousHelper<Eq> Setup;

BOOST_AUTO_TEST_CASE(BasicInsertionTest) {
	Setup::Vector test_r0({0.0, 0.1, 0.1, 1.0, 1.0});
	capd::ddeshelper::CubeSet<Setup::Vector> cubes(test_r0);
	cubes.insert({0,0,0});
	cubes.insert({0,1,0});
	cubes.insert({0,-1,0});
	cubes.insert({0,0,-1});
	cubes.insert({0,0,1});

    // Smoke test: iterate to ensure no crash
	for (auto it = cubes.csbegin(); it != cubes.csend(); ++it){
		// std::cout << (*it) << std::endl;
	}
	for (auto it = cubes.begin(); it != cubes.end(); ++it){
		// std::cout << (*it) << " " << Setup::Vector(it) << std::endl;
	}
    // No assertions here, just ensuring it runs as in the original test
}

BOOST_AUTO_TEST_CASE(InsertCoverGoodBoxTest) {
	Setup::Vector test_r0({0.0, 0.1, 0.1, 1.0, 1.0});
	Setup::Vector v(5);
	capd::ddeshelper::CubeSet<Setup::Vector> cubes2(test_r0);
	v[0] = Setup::Scalar(-0.2, 0.2);
	v[1] = Setup::Scalar(-0.3, 0.3);
	v[2] = Setup::Scalar(-0.2, 0.2);
	v[3] = Setup::Scalar(-0.03, 0.03);
	v[4] = Setup::Scalar(-0.5, 0.5);
	auto last_cut_needed = cubes2.insert_cover(v, 4);
	// cout << "last_cut_needed = " << last_cut_needed << ", should be smaller than 4" << std::endl;
    BOOST_CHECK_LT(last_cut_needed, 4);
}

BOOST_AUTO_TEST_CASE(InsertCoverTooBigBoxTest) {
	Setup::Vector test_r0({0.0, 0.1, 0.1, 1.0, 1.0});
	Setup::Vector v(5);
    // Initialize v properly
    v[0] = Setup::Scalar(-0.2, 0.2);
	v[1] = Setup::Scalar(-0.3, 0.3);
	v[2] = Setup::Scalar(-0.2, 0.2);
	v[3] = Setup::Scalar(-0.1, 5.1); // Modified from GoodBox
	v[4] = Setup::Scalar(-1.5, 0.5); // Modified from GoodBox

	capd::ddeshelper::CubeSet<Setup::Vector> cubes3(test_r0);
	auto last_cut_needed = cubes3.insert_cover(v, 4);
	// cout << "last_cut_needed = " << last_cut_needed << ", will be >= 4" << std::endl;
    BOOST_CHECK_GE(last_cut_needed, 4);

    // Reconstructing cubes2 for comparison (from GoodBoxTest)
    capd::ddeshelper::CubeSet<Setup::Vector> cubes2(test_r0);
    Setup::Vector v2(5);
	v2[0] = Setup::Scalar(-0.2, 0.2);
	v2[1] = Setup::Scalar(-0.3, 0.3);
	v2[2] = Setup::Scalar(-0.2, 0.2);
	v2[3] = Setup::Scalar(-0.03, 0.03);
	v2[4] = Setup::Scalar(-0.5, 0.5);
	cubes2.insert_cover(v2, 4);

    BOOST_CHECK(cubes3 == cubes3);
    BOOST_CHECK(cubes2 == cubes2);
    // They should be different
    BOOST_CHECK(!(cubes3 == cubes2));
}

BOOST_AUTO_TEST_CASE(StreamOperatorTest) {
	Setup::Vector test_r0({0.0, 0.1, 0.1, 1.0, 1.0});
    // Reconstruct 'cubes' from first test
	capd::ddeshelper::CubeSet<Setup::Vector> cubes(test_r0);
	cubes.insert({0,0,0});
	cubes.insert({0,1,0});
	cubes.insert({0,-1,0});
	cubes.insert({0,0,-1});
	cubes.insert({0,0,1});

	capd::ddeshelper::CubeSet<Setup::Vector> cubes_empty(test_r0);

    std::ostringstream oss;
    oss << cubes;
    BOOST_CHECK(!oss.str().empty());

    std::ostringstream oss_empty;
    oss_empty << cubes_empty;
    BOOST_CHECK(!oss_empty.str().empty());

	std::istringstream no_ws_input("{[0,-1,0],[0,0,-1],[0,0,0],[0,0,1],[0,1,0]} asadasd wsdaw wada");
	std::istringstream ws_input("{    \n   [0,-1,0]\n,[0,0,-1]\n     ,[0,0,0]     \n,[0,0,1]   \n    ,\t[0,1,0]  \t\n\n}\n\n wdadwa  w ada");
	capd::ddeshelper::CubeSet<Setup::Vector> cubes4(test_r0);
	no_ws_input >> cubes4;
	capd::ddeshelper::CubeSet<Setup::Vector> cubes5(test_r0);
	ws_input >> cubes5;

    BOOST_CHECK(cubes4 == cubes5);
    BOOST_CHECK(cubes == cubes4);
    BOOST_CHECK(cubes == cubes5);

	std::ostringstream osstest;
	osstest << cubes5 << "\n";
	std::istringstream isstest(osstest.str());
	capd::ddeshelper::CubeSet<Setup::Vector> cubes6(test_r0);
	isstest >> cubes6;
    BOOST_CHECK(cubes6 == cubes5);
}
