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
using TestGrid = DiscreteTimeGrid<RealType>;
using TestTimePoint = typename TestGrid::TimePointType;

BOOST_FIXTURE_TEST_CASE(GridConstructionTest, DiscreteTimeGridFixture<RealType>)
{
    BOOST_CHECK_EQUAL(grid.h(), step);
    
    // Test trivial grid
    TestGrid trivial;
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