#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE DDESolutionCurve
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

BOOST_AUTO_TEST_CASE( test_constructors ) {
    GridType grid(ScalarType(0.1));
    SolutionCurve curve(grid, 2, 2);
    BOOST_CHECK_EQUAL(curve.dimension(), 2);
    BOOST_CHECK_EQUAL(curve.storageN0(), 2);

    // Test t0, d, N0 constructor
    SolutionCurve curve2(grid.point(0), 3, 4);
    BOOST_CHECK_EQUAL(curve2.dimension(), 3);
    BOOST_CHECK_EQUAL(curve2.storageN0(), 4);
}

BOOST_AUTO_TEST_CASE( test_constructors_with_values ) {
    GridType grid(ScalarType(0.1));
    VectorType v(2); v[0] = 1.0; v[1] = 2.0;

    // Test t0, VectorType constructor
    SolutionCurve curve(grid.point(2), v);
    BOOST_CHECK_EQUAL(curve.dimension(), 2);
    BOOST_CHECK_EQUAL(curve.t0(), grid.point(2));
    BOOST_CHECK_EQUAL(curve.length(), 0);
    BOOST_CHECK_EQUAL(VectorType(curve.get_x()), v);

    // Test t0, SetType constructor
    SetType set_v(v);
    SolutionCurve curve2(grid.point(3), set_v);
    BOOST_CHECK_EQUAL(curve2.dimension(), 2);
    BOOST_CHECK_EQUAL(curve2.t0(), grid.point(3));
    BOOST_CHECK_EQUAL(VectorType(curve2.get_x()), v);
}

BOOST_AUTO_TEST_CASE( test_reduce_and_midcurve ) {
    GridType grid(ScalarType(0.1));
    VectorType v(2); v[0] = 1.0; v[1] = 2.0;

    SolutionCurve curve(grid.point(-2), grid.point(0), 2, v, 2);

    SolutionCurve mid = curve.midCurve();
    BOOST_CHECK_EQUAL(mid.length(), 2);
    BOOST_CHECK_EQUAL(mid.dimension(), 2);
}

BOOST_AUTO_TEST_CASE( test_reduce ) {
    GridType grid(ScalarType(0.1));
    VectorType v(2); v[0] = 1.0; v[1] = 2.0;

    SolutionCurve curve(grid.point(-2), grid.point(0), 3, v, 2);
    SolutionCurve reduced(grid.point(-2), grid.point(0), 2, v, 2);

    BOOST_CHECK_NO_THROW(curve.reduce(reduced));

    // Test copy and assignment
    SolutionCurve curve_copy(curve);
    BOOST_CHECK_EQUAL(curve_copy.length(), curve.length());

    SolutionCurve curve_assign = curve;
    BOOST_CHECK_EQUAL(curve_assign.length(), curve.length());
}

BOOST_AUTO_TEST_CASE( test_dot ) {
    GridType grid(ScalarType(0.1));
    VectorType v(2); v[0] = 1.0; v[1] = 2.0;
    SolutionCurve curve(grid.point(-2), grid.point(0), 3, v, 2);

    VectorType v2(2); v2[0] = 0.5; v2[1] = 1.0;

    // Test the dot product, returns double interval
    ScalarType result = curve.dot(v2);
    BOOST_CHECK_EQUAL(result, ScalarType(2.5)); // 1.0*0.5 + 2.0*1.0
}
