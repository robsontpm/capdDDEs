#define BOOST_TEST_MODULE BasicDiscreteDelaysFunctionalMapTests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <capd/ddes/BasicDiscreteDelaysFunctionalMap.hpp>
#include <capd/capdlib.h>
#include <vector>

// Mock types
using ScalarType = double;
using VectorType = capd::DVector;
using MatrixType = capd::DMatrix;
using TimePointType = double;

// Mock Jet
struct MockJet {
    std::vector<VectorType> coeffs;

    MockJet() {}
    explicit MockJet(const std::vector<VectorType>& c) : coeffs(c) {}

    int order() const { return coeffs.size() - 1; }
    int dimension() const { return coeffs.empty() ? 0 : coeffs[0].dimension(); }

    VectorType operator[](int k) const {
        if (k < 0 || k >= coeffs.size()) return VectorType(dimension());
        return coeffs[k];
    }

    VectorType& operator[](int k) {
        if (k >= coeffs.size()) {
            size_t new_size = k + 1;
            if (coeffs.empty()) {
                // Cannot resize without dimension, assume dimension 1 for empty mock if needed or handle logic
            } else {
               coeffs.resize(new_size, VectorType(dimension()));
            }
        }
        return coeffs[k];
    }
};

// Mock Solution Curve
struct MockSolutionCurve {    
    using VectorType = capd::DVector;
    using MatrixType = capd::DMatrix;
    using ScalarType = double;
    using TimePointType = double;
    using JetType = MockJet;
    using size_type = std::size_t;
    using DataType = VectorType;
    using RealType = ScalarType;

    int dim;

    explicit MockSolutionCurve(int d) : dim(d) {}

    int dimension() const { return dim; }

    VectorType eval(TimePointType t) const {
        VectorType v(dim);
        for(int i=0; i<dim; ++i) v[i] = t; // Simple function x(t) = t
        return v;
    }

    JetType jet(TimePointType t) const {
        // Return x(t) = t, x'(t) = 1, x''(t) = 0
        std::vector<VectorType> c(3);
        c[0] = eval(t);
        c[1] = VectorType(dim); for(int i=0; i<dim; ++i) c[1][i] = 1.0;
        c[2] = VectorType(dim); // 0
        return JetType(c);
    }

    VectorType getValueAtCurrent() const {
        return eval(0.0);
    }

    TimePointType rightDomain() const { return 0.0; }
};

// Mock Map
struct MockMap {
    int imgDim;
    int inputDim;

    explicit MockMap(int d, int inD) : imgDim(d), inputDim(inD) {}

    int imageDimension() const { return imgDim; }
    int dimension() const { return inputDim; }

    template<typename T, typename V>
    void operator()(const T& t, const V& in, V& out) const {
        // Simple map: out[i] = sum(in[j]) for all j
        for(int i=0; i<imgDim; ++i) {
            out[i] = T(0.0);
            for(int j=0; j<in.dimension(); ++j) {
                out[i] = out[i] + in[j];
            }
        }
    }
};

// Explicit instantiation to force compiler to generate code
// template class capd::ddes::BasicDiscreteDelaysFunctionalMap<MockMap, MockSolutionCurve, MockJet>;

// Test Suite
BOOST_AUTO_TEST_SUITE(BasicDiscreteDelaysFunctionalMapTests)

BOOST_AUTO_TEST_CASE(ConstructorTest) {
    int dim = 2;
    int delays = 1;
    MockMap map(dim, dim * (1 + delays));
    std::vector<double> delayVec = {1.0};

    capd::ddes::BasicDiscreteDelaysFunctionalMap<MockMap, MockSolutionCurve> funcMap(map, delayVec);

    BOOST_CHECK_EQUAL(funcMap.imageDimension(), dim);
    BOOST_CHECK_EQUAL(funcMap.dimension(), dim * (1 + delays));
    BOOST_CHECK_EQUAL(funcMap.delaysCount(), delays);
}

BOOST_AUTO_TEST_CASE(ConstructorNoDelaysTest) {
    int dim = 2;
    MockMap map(dim, dim);
    std::vector<double> delayVec = {};

    capd::ddes::BasicDiscreteDelaysFunctionalMap<MockMap, MockSolutionCurve> funcMap(map, delayVec);

    BOOST_CHECK_EQUAL(funcMap.delaysCount(), 0);
    BOOST_CHECK_EQUAL(funcMap.dimension(), dim);
}

BOOST_AUTO_TEST_CASE(DimensionMismatchTest) {
    int dim = 2;
    MockMap map(dim, 10);
    std::vector<double> delayVec = {1.0};

    BOOST_CHECK_THROW(
        (capd::ddes::BasicDiscreteDelaysFunctionalMap<MockMap, MockSolutionCurve>(map, delayVec)),
        std::logic_error
    );
}

BOOST_AUTO_TEST_CASE(EvaluationTest) {
    int dim = 1;
    int delays = 1;
    // Map: x(t) + x(t-1)
    // x(t) = t.
    // t=2. x(2)=2, x(1)=1. Sum = 3.
    MockMap map(dim, dim * (1 + delays));
    std::vector<double> delayVec = {1.0};

    capd::ddes::BasicDiscreteDelaysFunctionalMap<MockMap, MockSolutionCurve> funcMap(map, delayVec);
    MockSolutionCurve curve(dim);

    capd::DVector result = funcMap(2.0, curve);

    BOOST_CHECK_EQUAL(result[0], 3.0);
}

BOOST_AUTO_TEST_CASE(CollectComputationDataTest) {
    int dim = 1;
    int delays = 1;
    MockMap map(dim, dim * (1 + delays));
    std::vector<double> delayVec = {1.0};

    capd::ddes::BasicDiscreteDelaysFunctionalMap<MockMap, MockSolutionCurve> funcMap(map, delayVec);
    MockSolutionCurve curve(dim);

    capd::ddes::BasicDiscreteDelaysFunctionalMap<MockMap, MockSolutionCurve>::VariableStorageType out_u;
    capd::ddes::BasicDiscreteDelaysFunctionalMap<MockMap, MockSolutionCurve>::size_type order;

    funcMap.collectComputationData(0.0, 0.1, 0.1, curve, out_u, order);

    // out_u should contain: x(0), x(-1), x'(-1), x''(-1) ... depending on order
    // x(-1) = -1. x'(-1) = 1.
    // CollectComputationData:
    // jets.push_back(x.jet(t0 - *tau)); -> jet at -1. order=2.
    // out_admissible_order becomes 3.
    // out_u:
    // 1. x.getValueAtCurrent() -> x(0) = 0.
    // 2. jets loop:
    // k=0: jet at -1 [0] -> -1
    // k=1: jet at -1 [1] -> 1
    // k=2: jet at -1 [2] -> 0

    // Total size should be 1 + delays * admissible_order = 1 + 1 * 3 = 4.

    BOOST_CHECK_EQUAL(out_u.size(), 4);
    BOOST_CHECK_EQUAL(out_u[0][0], 0.0); // x(0)
    BOOST_CHECK_EQUAL(out_u[1][0], -1.0); // x(-1)
    BOOST_CHECK_EQUAL(out_u[2][0], 1.0); // x'(-1)
    BOOST_CHECK_EQUAL(out_u[3][0], 0.0); // x''(-1)
}

BOOST_AUTO_TEST_CASE(MakeCompatibleSegmentTest) {
    int dim = 1;
    MockMap map(dim, dim);
    std::vector<double> delayVec;

    capd::ddes::BasicDiscreteDelaysFunctionalMap<MockMap, MockSolutionCurve> funcMap(map, delayVec);
    capd::DVector v(dim);

    BOOST_CHECK_THROW(funcMap.makeCompatibleSegment(v), std::logic_error);
}

BOOST_AUTO_TEST_CASE(ComputeDDECoefficientsTest) {
    int dim = 1;
    int delays = 1;
    MockMap map(dim, dim * (1 + delays));
    std::vector<double> delayVec = {1.0};

    capd::ddes::BasicDiscreteDelaysFunctionalMap<MockMap, MockSolutionCurve> funcMap(map, delayVec);

    // Setup u
    // u[0] = x(t0)
    // u[1] = x(t0-tau)
    capd::ddes::BasicDiscreteDelaysFunctionalMap<MockMap, MockSolutionCurve>::ValueStorageType u(2);
    u[0] = capd::DVector(dim); u[0][0] = 1.0;
    u[1] = capd::DVector(dim); u[1][0] = 2.0;

    capd::ddes::BasicDiscreteDelaysFunctionalMap<MockMap, MockSolutionCurve>::ValueStorageType coeffs(2); // Order 1 (coeffs[0], coeffs[1])

    funcMap.computeDDECoefficients(0.0, u, coeffs);

    // Check results
    // coeffs[0] = u[0] = 1.0
    // coeffs[1] = f(u) / 1
    // f(u) = u[0] + u[1] = 1 + 2 = 3.
    // So coeffs[1] should be 3.

    BOOST_CHECK_EQUAL(coeffs[0][0], 1.0);
    BOOST_CHECK_EQUAL(coeffs[1][0], 3.0);
}

BOOST_AUTO_TEST_CASE(ComputeDDECoefficientsWithJacobianTest) {
    int dim = 1;
    int delays = 1;
    MockMap map(dim, dim * (1 + delays));
    std::vector<double> delayVec = {1.0};

    capd::ddes::BasicDiscreteDelaysFunctionalMap<MockMap, MockSolutionCurve> funcMap(map, delayVec);

    capd::ddes::BasicDiscreteDelaysFunctionalMap<MockMap, MockSolutionCurve>::ValueStorageType u(2);
    u[0] = capd::DVector(dim); u[0][0] = 1.0;
    u[1] = capd::DVector(dim); u[1][0] = 2.0;

    capd::ddes::BasicDiscreteDelaysFunctionalMap<MockMap, MockSolutionCurve>::ValueStorageType coeffs(2);
    capd::ddes::BasicDiscreteDelaysFunctionalMap<MockMap, MockSolutionCurve>::JacobianStorageType Du;

    funcMap.computeDDECoefficients(0.0, u, coeffs, Du);

    // Check coeffs
    BOOST_CHECK_EQUAL(coeffs[0][0], 1.0);
    BOOST_CHECK_EQUAL(coeffs[1][0], 3.0);

    // Check Jacobian Du
    // Du[k][j] = d(coeffs[k]) / d(u[j])
    // k=0: coeffs[0] = u[0].
    // d(coeffs[0])/d(u[0]) = Id = 1.
    // d(coeffs[0])/d(u[1]) = 0.

    BOOST_CHECK_EQUAL(Du[0][0][0][0], 1.0);
    BOOST_CHECK_EQUAL(Du[0][1][0][0], 0.0);

    // k=1: coeffs[1] = u[0] + u[1]
    // d(coeffs[1])/d(u[0]) = 1.
    // d(coeffs[1])/d(u[1]) = 1.

    BOOST_CHECK_EQUAL(Du[1][0][0][0], 1.0);
    BOOST_CHECK_EQUAL(Du[1][1][0][0], 1.0);
}

BOOST_AUTO_TEST_SUITE_END()
