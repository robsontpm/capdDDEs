#define BOOST_TEST_MODULE DDECommonTests
#include <boost/test/included/unit_test.hpp>
#include <capd/ddes/DDECommon.h>
#include <capd/capdlib.h>

using namespace capd::ddes;

// Test fixture for DiscreteTimeGrid tests
template<typename RealType>
struct DiscreteTimeGridFixture {
    const RealType step;    
    DiscreteTimeGrid<RealType> grid;

    DiscreteTimeGridFixture() : step(1.0), grid(step) {}
};

BOOST_AUTO_TEST_SUITE(HelperFunctionTests)

BOOST_AUTO_TEST_CASE(SafeDeleteTest)
{
    int* ptr = new int(42);
    bool is_owner = true;
    
    helper_safe_delete(ptr, is_owner);
    BOOST_CHECK_EQUAL(ptr, nullptr);

    // Test non-owner case
    int* ptr2 = new int(42);
    is_owner = false;
    helper_safe_delete(ptr2, is_owner);
    BOOST_CHECK_NE(ptr2, nullptr);
    delete ptr2; // Clean up
}

BOOST_AUTO_TEST_CASE(SafeArrayDeleteTest)
{
    int* arr = new int[5];
    bool is_owner = true;
    
    helper_safe_array_delete(arr, is_owner);
    BOOST_CHECK_EQUAL(arr, nullptr);

    // Test non-owner case
    int* arr2 = new int[5];
    is_owner = false;
    helper_safe_array_delete(arr2, is_owner);
    BOOST_CHECK_NE(arr2, nullptr);
    delete[] arr2; // Clean up
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(DiscreteTimeGridTests)

using RealType = double; // You can change this to test with different types
using Grid = DiscreteTimeGrid<RealType>;
using TimePoint = typename Grid::TimePointType;

BOOST_FIXTURE_TEST_CASE(GridConstructionTest, DiscreteTimeGridFixture<RealType>)
{
    BOOST_CHECK_EQUAL(grid.h(), step);
    
    // Test trivial grid
    Grid trivial;
    BOOST_CHECK_EQUAL(trivial.h(), 0.0);
}

BOOST_FIXTURE_TEST_CASE(TimePointCreationTest, DiscreteTimeGridFixture<RealType>)
{
    auto point = grid.point(5);
    BOOST_CHECK_EQUAL(static_cast<RealType>(point), 5.0 * step);
    BOOST_CHECK_EQUAL(point.toInt(), 5);
}

BOOST_FIXTURE_TEST_CASE(TimePointArithmeticTest, DiscreteTimeGridFixture<RealType>)
{
    auto p1 = grid.point(5);
    auto p2 = grid.point(3);
    
    // Addition
    auto sum = p1 + p2;
    BOOST_CHECK_EQUAL(sum.toInt(), 8);
    
    // Subtraction
    auto diff = p1 - p2;
    BOOST_CHECK_EQUAL(diff.toInt(), 2);
    
    // Increment/Decrement
    ++p1;
    BOOST_CHECK_EQUAL(p1.toInt(), 6);
    --p1;
    BOOST_CHECK_EQUAL(p1.toInt(), 5);
}

BOOST_FIXTURE_TEST_CASE(TimePointComparisonTest, DiscreteTimeGridFixture<RealType>)
{
    auto p1 = grid.point(5);
    auto p2 = grid.point(3);
    
    BOOST_CHECK(p1 > p2);
    BOOST_CHECK(p2 < p1);
    BOOST_CHECK(p1 >= p2);
    BOOST_CHECK(p2 <= p1);
    BOOST_CHECK(p1 != p2);
    
    auto p3 = grid.point(5);
    BOOST_CHECK(p1 == p3);
}

// Test the default constructor
BOOST_AUTO_TEST_CASE(DefaultConstructorTest) {
    Grid grid;
    BOOST_CHECK_EQUAL(grid.h(), 0);
}

// Test the constructor with step size
BOOST_AUTO_TEST_CASE(ConstructorWithStepSizeTest) {
    double h = 0.1;
    Grid grid(h);
    BOOST_CHECK_EQUAL(grid.h(), h);
}

// Test the point creation
BOOST_AUTO_TEST_CASE(PointCreationTest) {
    double h = 0.1;
    Grid grid(h);
    TimePoint point = grid.point(5);
    BOOST_CHECK_EQUAL(static_cast<double>(point), h * 5);
}

// Test the operator==
BOOST_AUTO_TEST_CASE(EqualityOperatorTest) {
    double h = 0.1;
    Grid grid1(h);
    Grid grid2(grid1);
    // only the copy of a given grid works! See the docs for default constructor!
    BOOST_CHECK(grid1 == grid2);
}

// Test the operator!=
BOOST_AUTO_TEST_CASE(InequalityOperatorTest) {
    double h1 = 0.1;
    double h2 = 0.2;
    Grid grid1(h1);
    Grid grid2(h2);
    BOOST_CHECK(grid1 != grid2);
}

// Test the TimePoint addition
BOOST_AUTO_TEST_CASE(TimePointAdditionTest) {
    double h = 0.1;
    Grid grid(h);
    TimePoint point1 = grid.point(5);
    TimePoint point2 = grid.point(3);
    TimePoint result = point1 + point2;
    BOOST_CHECK_EQUAL(static_cast<double>(result), h * 8);
}

// Test the TimePoint subtraction
BOOST_AUTO_TEST_CASE(TimePointSubtractionTest) {
    double h = 0.1;
    Grid grid(h);
    TimePoint point1 = grid.point(5);
    TimePoint point2 = grid.point(3);
    TimePoint result = point1 - point2;
    BOOST_CHECK_EQUAL(static_cast<double>(result), h * 2);
}

// Test the TimePoint increment
BOOST_AUTO_TEST_CASE(TimePointIncrementTest) {
    double h = 0.1;
    Grid grid(h);
    TimePoint point = grid.point(5);
    ++point;
    BOOST_CHECK_EQUAL(static_cast<double>(point), h * 6);
}

// Test the TimePoint decrement
BOOST_AUTO_TEST_CASE(TimePointDecrementTest) {
    double h = 0.1;
    Grid grid(h);
    TimePoint point = grid.point(5);
    --point;
    BOOST_CHECK_EQUAL(static_cast<double>(point), h * 4);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(TaylorSumTests)

BOOST_AUTO_TEST_CASE(SumTaylorForwardTest)
{
    std::vector<double> coefficients = {1.0, 2.0, 3.0}; // represents 1 + 2x + 3 x^2
    double step = 2.0;
    double result = 0.0;
    
    sumTaylorForward(coefficients.begin(), 2, step, result);
    // Expected: 1 + 2(2) + 3 (4) = 17
    BOOST_CHECK_CLOSE(result, 17.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(SumTaylorBackwardTest)
{
    std::vector<double> coefficients = {1.0, 2.0, 3.0}; // represents 1 + 2x + 3x^2
    double step = 2.0;
    double result = 0.0;
    
    sumTaylorBackward(coefficients.end() - 1, 2, step, result);
    // Expected: 1 + 2(2) + 3(4) = 17
    BOOST_CHECK_CLOSE(result, 17.0, 1e-10);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(MatrixOperationsTests)

BOOST_AUTO_TEST_CASE(ExtractDiagonalBlocksTest)
{
    // Create a test matrix
    capd::DMatrix M(4, 4);
    // Fill with test data
    for(size_t i = 0; i < 4; ++i)
        for(size_t j = 0; j < 4; ++j)
            M[i][j] = i * 4 + j;
            
    size_t offBlockCount = 0;
    auto blocks = extractDiagonalBlocks(M, size_t(2), offBlockCount);
    
    BOOST_CHECK_EQUAL(blocks.size(), 2);
    BOOST_CHECK_EQUAL(blocks[0].numberOfRows(), 2);
    BOOST_CHECK_EQUAL(blocks[0].numberOfColumns(), 2);
    
    // Test invalid block size
    BOOST_CHECK_THROW(extractDiagonalBlocks(M, size_t(3), offBlockCount), std::range_error);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(DDECommonAdditionalTests)

BOOST_AUTO_TEST_CASE(HelperDumpTests)
{
    std::stringstream ss;
    ss << "line1\nline2";
    helper_dump_line(ss);
    // Should have read "line1"
    std::string remaining;
    ss >> remaining;
    BOOST_CHECK_EQUAL(remaining, "line2");

    std::stringstream ss2;
    ss2 << "badge1 badge2";
    helper_dump_badge(ss2);
    // Should have read "badge1"
    ss2 >> remaining;
    BOOST_CHECK_EQUAL(remaining, "badge2");
}

BOOST_AUTO_TEST_CASE(EcloseStepTest)
{
    // Test double version
    double h = 0.5;
    BOOST_CHECK_EQUAL(ecloseStep(h), 0.5);

    // Test Interval version
    capd::interval h_int(0.5);
    capd::interval res = ecloseStep(h_int);
    // res should be [0, 1] * 0.5 = [0, 0.5]
    BOOST_CHECK_EQUAL(res.leftBound(), 0.0);
    BOOST_CHECK_EQUAL(res.rightBound(), 0.5);
}

BOOST_AUTO_TEST_CASE(ShowEnclosedIntervalTest)
{
    // Test double version
    std::string s = showEnclosedInterval(1.0, 2.0);
    BOOST_CHECK_EQUAL(s, "[1, 2)");

    // Test Interval version
    capd::interval a(1.0);
    capd::interval b(2.0);
    std::string s2 = showEnclosedInterval(a, b);
    BOOST_CHECK_EQUAL(s2, "[1, 2)");
}

BOOST_AUTO_TEST_CASE(RethrowTest)
{
    try {
        throw std::runtime_error("original error");
    } catch (const std::exception& e) {
        // Test rethrow<std::exception>
        // rethrow returns the exception object, it doesn't throw it.
        std::runtime_error new_e = rethrow<std::exception, std::runtime_error>("Context", e);
        // glue is "\n    "
        BOOST_CHECK(std::string(new_e.what()).find("Context\n    original error") != std::string::npos);

        // Test rethrow<std::exception, std::logic_error>
        std::logic_error new_logic = rethrow<std::exception, std::logic_error>("Logic", e);
        BOOST_CHECK(std::string(new_logic.what()).find("Logic\n    original error") != std::string::npos);
    }
}

BOOST_AUTO_TEST_CASE(ClosestIntTest)
{
    // Test generic version (double)
    BOOST_CHECK_EQUAL(closestInt(3.1), 3);
    BOOST_CHECK_EQUAL(closestInt(3.9), 3); // int(3.9) is 3
    BOOST_CHECK_EQUAL(closestInt(-3.1), -3);
    BOOST_CHECK_EQUAL(closestInt(-3.9), -3);

    // Test closestSmallerInt
    BOOST_CHECK_EQUAL(closestSmallerInt(3.1), 3);
    BOOST_CHECK_EQUAL(closestSmallerInt(3.9), 3);
    BOOST_CHECK_EQUAL(closestSmallerInt(-3.1), -4);

    // Test capd::interval version
    capd::interval iv(3.5);
    BOOST_CHECK_EQUAL(closestInt(iv), 3);
    BOOST_CHECK_EQUAL(closestSmallerInt(iv), 3);

    capd::interval iv_neg(-3.5);
    BOOST_CHECK_EQUAL(closestInt(iv_neg), -3);
    BOOST_CHECK_EQUAL(closestSmallerInt(iv_neg), -4);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(DiscreteTimeGridExtendedTests)

using RealType = double;
using Grid = DiscreteTimeGrid<RealType>;
using TimePoint = typename Grid::TimePointType;

BOOST_AUTO_TEST_CASE(SplitTest)
{
    double h = 0.5;
    Grid grid(h);

    // Test exact point
    double t = 1.0; // 2 * h
    TimePoint ti = grid.point(0);
    double epsi;
    grid.split(t, ti, epsi);

    BOOST_CHECK_EQUAL(ti.toInt(), 2);
    BOOST_CHECK_SMALL(epsi, 1e-14);

    // Test point slightly after grid point
    t = 1.1; // 2*h + 0.1
    grid.split(t, ti, epsi);
    BOOST_CHECK_EQUAL(ti.toInt(), 2);
    BOOST_CHECK_CLOSE(epsi, 0.1, 1e-10);

    // Test point slightly before grid point
    t = 0.9; // 1.8*h
    grid.split(t, ti, epsi);
    BOOST_CHECK_EQUAL(ti.toInt(), 1);
    BOOST_CHECK_CLOSE(epsi, 0.4, 1e-10);

    // Let's test with negative numbers
    t = -0.9;
    grid.split(t, ti, epsi);
    BOOST_CHECK_EQUAL(ti.toInt(), -2);
    BOOST_CHECK_CLOSE(epsi, 0.1, 1e-10);
}

BOOST_AUTO_TEST_CASE(TimePointExtendedTest)
{
    Grid grid(0.1);
    TimePoint p = grid.point(0);

    // Test isZero
    BOOST_CHECK(p.isZero());
    TimePoint p2 = grid.point(1);
    BOOST_CHECK(!p2.isZero());

    // Test show
    // format: value := badge i h
    // 0.1 := DiscreteTimePoint 1 0.1
    std::string s = p2.show();
    BOOST_CHECK(s.find("0.1") != std::string::npos);
    BOOST_CHECK(s.find("1") != std::string::npos);

    // Test stream IO
    std::stringstream ss;
    ss << p2;
    // Format: value := badge i h
    // We must initialize p_target with compatible grid
    TimePoint p_target = grid.point(0);
    ss >> p_target;
    BOOST_CHECK(p_target == p2);
}

BOOST_AUTO_TEST_CASE(GridCompatibilityTest)
{
    Grid g1(0.1);
    Grid g2(0.2); // Different h
    Grid g3(0.1); // Same h, but different object.

    TimePoint p1 = g1.point(1);
    TimePoint p2 = g2.point(1);
    TimePoint p3 = g3.point(1);

    // Arithmetic with incompatible grids should throw
    BOOST_CHECK_THROW(p1 + p2, std::logic_error);
    BOOST_CHECK_THROW(p1 + p3, std::logic_error);

    // Comparison
    // p1 < p2 falls back to value comparison
    BOOST_CHECK(p1 < p2);

    // Zero compatibility
    TimePoint z1 = g1.point(0);
    TimePoint z2 = g2.point(0);

    TimePoint sum = p1 + z2;
    BOOST_CHECK_EQUAL(sum.toInt(), 1);
    BOOST_CHECK(sum.sameGrid(p1));
}

BOOST_AUTO_TEST_CASE(ExtractDiagonalBlocksExtendedTest)
{
    capd::DMatrix M(5, 5); // 5x5
    size_t off;
    // extract with block size 2.
    // 5 is not multiple of 2. Should throw range_error.
    BOOST_CHECK_THROW(extractDiagonalBlocks(M, size_t(2), off), std::range_error);

    capd::DMatrix M4(4, 4);
    // Fill diagonal with blocks
    // Block 0: rows 0-1, cols 0-1.
    M4[0][0] = 1; M4[1][1] = 1;
    // Block 1: rows 2-3, cols 2-3.
    M4[2][2] = 1; M4[3][3] = 1;
    // Off diagonal element
    M4[0][3] = 9;

    off = 0;
    auto blocks = extractDiagonalBlocks(M4, size_t(2), off);
    BOOST_CHECK_EQUAL(blocks.size(), 2);
    // off should be 1
    BOOST_CHECK_EQUAL(off, 1);
}

BOOST_AUTO_TEST_SUITE_END()
