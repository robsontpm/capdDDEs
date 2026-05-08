#define BOOST_TEST_MODULE DDESolutionCurveBugTests3
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <capd/ddes/DDESolutionCurve.h>
#include <capd/ddes/DDESolutionCurve.hpp>
#include <capd/ddes/DDECommon.h>
#include <capd/ddes/storage/SharedDoubleton.h>
#include <capd/ddes/storage/GenericJet.h>
#include <capd/ddes/storage/GenericJet.hpp>
#include <capd/capdlib.h>

using namespace capd::ddes;

typedef capd::DInterval Interval;
typedef capd::IMatrix IMatrix;
typedef capd::ddes::SharedDoubleton<IMatrix> DDESet;
typedef DiscreteTimeGrid<Interval> GridType;
typedef GridType::TimePointType TimePointType;

BOOST_AUTO_TEST_SUITE(DDESolutionCurve_BUG_Suite)

BOOST_AUTO_TEST_CASE(Reduce_Exception) {
    GridType grid(Interval(0.1));
    TimePointType t0 = grid.point(0);
    TimePointType t1 = grid.point(2);

    capd::IVector vec1(2);
    vec1[0]=3.0; vec1[1]=4.0;

    DDESolutionCurve<DDESet> curve_high(t0, t1, 3, vec1, 3);
    DDESolutionCurve<DDESet> curve_low(t0, t1, 1, vec1, 3);

    BOOST_CHECK_THROW(curve_low.reduce(curve_high), std::logic_error);
}

BOOST_AUTO_TEST_SUITE_END()
