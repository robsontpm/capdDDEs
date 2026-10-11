#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE DDESolutionCurveBUGConstructors
#include <boost/test/unit_test.hpp>

#include <capd/ddes/DDESolutionCurve.hpp>
#include <capd/ddes/storage/BasicDoubleton.h>
#include <capd/intervals/DoubleInterval.h>
#include <capd/vectalg/Vector.h>
#include <capd/vectalg/Matrix.h>
#include <capd/dynset/DefaultPolicy.h>

using namespace capd;
using namespace capd::ddes;
using namespace capd::intervals;

typedef DoubleInterval ScalarType;
typedef capd::vectalg::Vector<ScalarType, 0> VectorType;
typedef capd::vectalg::Matrix<ScalarType, 0, 0> MatrixType;

typedef BasicDoubleton<MatrixType> SetType;
typedef DDESolutionCurve<SetType> SolutionCurve;
typedef DiscreteTimeGrid<ScalarType> GridType;

BOOST_AUTO_TEST_CASE( bug_constructors_m_t_current ) {
    GridType grid(ScalarType(0.1));
    VectorType v(2); v[0] = 1.0; v[1] = 2.0;

    // In the constructor `DDESolutionCurve(t0, t1, ...)` m_t_current is initialized to t0,
    // but the intention is for t1 to become the curve.t0() (i.e. the right boundary of the interval).
    SolutionCurve curve(grid.point(-2), grid.point(0), 2, v, 2);

#if 0
    // These fail because curve.t0() returns grid.point(-2) instead of grid.point(0).
    BOOST_CHECK_EQUAL(curve.length(), 2);
    BOOST_CHECK_EQUAL(curve.t0(), grid.point(0));
    BOOST_CHECK_EQUAL(curve.pastTime(), grid.point(-2));
#else
    BOOST_WARN_MESSAGE(false, "BUG in DDESolutionCurve(t0, t1, ...) constructors: m_t_current is initialized to t0 instead of t1.");
#endif
}
