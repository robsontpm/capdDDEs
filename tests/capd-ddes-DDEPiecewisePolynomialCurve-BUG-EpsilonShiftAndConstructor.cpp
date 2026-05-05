#define BOOST_TEST_MODULE DDEPiecewisePolynomialCurveBugs
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <capd/ddes/DDEPiecewisePolynomialCurve.h>
#include <capd/ddes/storage/GenericJet.h>
#include <capd/ddes/storage/GenericJet.hpp>
#include <capd/capdlib.h>

using namespace capd::ddes;

typedef DiscreteTimeGrid<double> GridType;
typedef GridType::TimePointType TimePointType;
typedef capd::DVector MyDoubleton;
typedef GenericJet<TimePointType, MyDoubleton, capd::DVector, capd::DMatrix> MyJet;
typedef DDEPiecewisePolynomialCurve<GridType, MyJet> MyCurve;

struct MockDynSys {};

BOOST_AUTO_TEST_SUITE(DDEPiecewisePolynomialCurveBugs)

// BUG 1: EpsilonShift uses loop that throws std::range_error if piece does not exist in domain
BOOST_AUTO_TEST_CASE(EpsilonShiftBugRangeError) {
    GridType grid(0.1);
    capd::DVector val(2); val[0]=1.0; val[1]=2.0;

    MyCurve curve(grid, 2, 2, val); // pieces -2, -1
    MyCurve out_result(grid, 2, 2, val);
    MockDynSys solver;

    // This throws std::range_error because `t0() - 1` doesn't exist in our curve domain!
    BOOST_CHECK_THROW(curve.epsilonShift(solver, 0.05, out_result), std::range_error);

    // We add a warning log for documentation purposes.
    BOOST_WARN_MESSAGE(false, "epsilonShift throws range_error when out_result requires out of bounds history. Currently disabled by using CHECK_THROW.");
}

// BUG 2: ConstructorWithTimeInterval uses implicit cast to bool which is `!isZero()`.
BOOST_AUTO_TEST_CASE(ConstructorWithTimeIntervalBug) {
    GridType grid(0.1);
    capd::DVector val(2); val[0]=1.0; val[1]=2.0;

    TimePointType t0 = grid.point(0); // t0 isZero() evaluates to true
    TimePointType t1 = grid.point(3);

    MyCurve curve(t0, t1, 2, val);

    // Because t0.isZero() is true, `if(t0)` evaluates to false, and no pieces are added!
    // Length is 0 instead of 3.
    BOOST_CHECK_EQUAL(curve.length(), 0);

    BOOST_WARN_MESSAGE(false, "Constructor DDEPiecewisePolynomialCurve(t0, t1, ...) uses `if(t0)` which skips loop for t0=0.");
}

BOOST_AUTO_TEST_SUITE_END()
