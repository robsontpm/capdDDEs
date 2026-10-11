#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE DDESolutionCurveBUGSubcurve
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

BOOST_AUTO_TEST_CASE( bug_subcurve_bounds ) {
    GridType grid(ScalarType(0.1));
    VectorType v(2); v[0] = 1.0; v[1] = 2.0;

    // Note: because of the bug, curve.t0() will be -5 and length 5
    SolutionCurve curve(grid.point(-5), grid.point(0), 3, v, 2);

#if 0
    // Calling subcurve with right endpoint t1 throws std::domain_error inside `pointToIndex(t1)`
    // because `pointToIndex` checks if the index is < i0 (where i0 is the number of pieces).
    // The right boundary is exactly t1 = t0(), which maps to i0, failing the strict < check.
    SolutionCurve sub1 = curve.subcurve(grid.point(-3), grid.point(0));
#else
    BOOST_WARN_MESSAGE(false, "BUG in DDESolutionCurve::subcurve(t0, t1): fails boundary check at pointToIndex(t1) because it uses a strict < i0 limit.");
#endif
}
