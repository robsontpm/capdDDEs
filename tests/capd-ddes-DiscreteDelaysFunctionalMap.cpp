#define BOOST_TEST_MODULE DiscreteDelaysFunctionalMapTests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <capd/ddes/DiscreteDelaysFunctionalMap.h>
#include <capd/ddes/DiscreteDelaysFunctionalMap.hpp>
#include <capd/capdlib.h>
#include <vector>

// Use Interval types for Rigorous Map
using ScalarType = capd::Interval;
using VectorType = capd::IVector;
using MatrixType = capd::IMatrix;
using TimePointType = capd::Interval;

// Mock Jet for Interval
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
                // Cannot resize without dimension, assume dimension 1 for empty mock if needed
            } else {
               coeffs.resize(new_size, VectorType(dimension()));
            }
        }
        return coeffs[k];
    }

    // Evaluate Taylor sum at delta: sum(c[k] * delta^k)
    VectorType evalAtDelta(const ScalarType& delta) const {
        VectorType result(dimension()); // Initialize with 0
        if (coeffs.empty()) return result;

        // Horner's method would be better but forward is fine for test
        // sum_{k=0}^n c_k * delta^k
        ScalarType power = 1.0;
        for(const auto& c : coeffs) {
            result = result + c * power;
            power = power * delta;
        }
        return result;
    }

    // Evaluate k-th derivative at delta: sum_{j=k}^n c_j * (j!/(j-k)!) * delta^(j-k)
    VectorType evalCoeffAtDelta(int k, const ScalarType& delta) const {
        VectorType result(dimension());
        if (k > order()) return result;

        ScalarType power = 1.0;
        for (int j = k; j <= order(); ++j) {
            // Binomial coeff (j choose k)
            double binom = 1.0;
            for(int i=0; i<k; ++i) binom = binom * (j - i) / (i + 1);
            // Factorial k!
            double k_fact = 1.0;
            for(int i=1; i<=k; ++i) k_fact *= i;

            // The formula: x^(k)(t) = sum_{j=k}^n c_j * (j!/(j-k)!) * delta^(j-k)
            // But evalCoeffAtDelta usually returns the k-th coefficient of the expansion at t0+delta, which is x^(k)(t0+delta)/k!.
            // So we need (1/k!) * sum_{j=k}^n c_j * (j!/(j-k)!) * delta^(j-k)
            // = sum_{j=k}^n c_j * binom(j, k) * delta^(j-k)

            result = result + coeffs[j] * (binom * power);
            power = power * delta;
        }
        return result;
    }
};

// Mock Solution Curve
struct MockSolutionCurve {
    using VectorType = capd::IVector;
    using MatrixType = capd::IMatrix;
    using ScalarType = capd::Interval;
    using TimePointType = capd::Interval;
    using JetType = MockJet;
    // Use unsigned int instead of std::size_t to match IVector::dimension() return type
    // and avoid template ambiguity in checkDimension.
    using size_type = unsigned int;
    using DataType = VectorType;
    using RealType = ScalarType;

    int dim;

    explicit MockSolutionCurve(int d) : dim(d) {}

    int dimension() const { return dim; }

    VectorType eval(TimePointType t) const {
        VectorType v(dim);
        // x(t) = t
        for(int i=0; i<dim; ++i) v[i] = t;
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

// Mock Map compatible with fadbad
struct MockMap {
    int imgDim;
    int inputDim;

    explicit MockMap(int d, int inD) : imgDim(d), inputDim(inD) {}

    int imageDimension() const { return imgDim; }
    int dimension() const { return inputDim; }

    template<typename T, typename V>
    void operator()(const T& t, const V& in, V& out) const {
        // Simple map: out[i] = sum(in[j]) for all j
        // Works for Interval, fadbad::T<Interval>, etc.
        for(int i=0; i<imgDim; ++i) {
            out[i] = T(0.0);
            for(int j=0; j<in.dimension(); ++j) {
                out[i] = out[i] + in[j];
            }
        }
    }
};

BOOST_AUTO_TEST_SUITE(DiscreteDelaysFunctionalMapTests)

BOOST_AUTO_TEST_CASE(ConstructorTest) {
    int dim = 2;
    int delays = 1;
    MockMap map(dim, dim * (1 + delays));
    std::vector<capd::Interval> delayVec = {1.0};

    capd::ddes::DiscreteDelaysFunctionalMap<MockMap, MockSolutionCurve> funcMap(map, delayVec);

    BOOST_CHECK_EQUAL(funcMap.imageDimension(), dim);
    BOOST_CHECK_EQUAL(funcMap.dimension(), dim * (1 + delays));
    BOOST_CHECK_EQUAL(funcMap.delaysCount(), delays);
}

BOOST_AUTO_TEST_CASE(EvaluationTest) {
    int dim = 1;
    int delays = 1;
    // Map: x(t) + x(t-1)
    // x(t) = t.
    // t=2. x(2)=2, x(1)=1. Sum = 3.
    MockMap map(dim, dim * (1 + delays));
    std::vector<capd::Interval> delayVec = {1.0};

    capd::ddes::DiscreteDelaysFunctionalMap<MockMap, MockSolutionCurve> funcMap(map, delayVec);
    MockSolutionCurve curve(dim);

    capd::IVector result = funcMap(capd::Interval(2.0), curve);

    BOOST_CHECK_CLOSE(result[0].mid(), 3.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(CollectComputationDataTest) {
    int dim = 1;
    int delays = 1;
    MockMap map(dim, dim * (1 + delays));
    std::vector<capd::Interval> delayVec = {1.0};

    capd::ddes::DiscreteDelaysFunctionalMap<MockMap, MockSolutionCurve> funcMap(map, delayVec);
    MockSolutionCurve curve(dim);

    capd::ddes::DiscreteDelaysFunctionalMap<MockMap, MockSolutionCurve>::VariableStorageType out_u;
    capd::ddes::DiscreteDelaysFunctionalMap<MockMap, MockSolutionCurve>::ValueStorageType out_encl;
    capd::ddes::DiscreteDelaysFunctionalMap<MockMap, MockSolutionCurve>::size_type order;

    // t0=0, th=0.1, dt=0.1
    // We use small dt so findRoughEnclosure can succeed easily.
    funcMap.collectComputationData(0.0, 0.1, 0.1, curve, out_u, out_encl, order);

    // out_u should contain: x(0), x(-1), x'(-1), x''(-1) ... depending on order
    // Order is determined by jets. back().order() + 1.
    // MockJet returns order 2 (size 3).
    // So order should be 3.
    // out_u size: 1 (x(0)) + 3 (jet components) = 4.

    BOOST_CHECK_EQUAL(out_u.size(), 4);
    BOOST_CHECK_CLOSE(out_u[0][0].mid(), 0.0, 1e-10); // x(0)
    BOOST_CHECK_CLOSE(out_u[1][0].mid(), -1.0, 1e-10); // x(-1)
}

BOOST_AUTO_TEST_CASE(ComputeDDECoefficientsTest) {
    int dim = 1;
    int delays = 1;
    MockMap map(dim, dim * (1 + delays));
    std::vector<capd::Interval> delayVec = {1.0};

    capd::ddes::DiscreteDelaysFunctionalMap<MockMap, MockSolutionCurve> funcMap(map, delayVec);

    // Setup u
    // u[0] = x(t0)
    // u[1] = x(t0-tau)
    capd::ddes::DiscreteDelaysFunctionalMap<MockMap, MockSolutionCurve>::ValueStorageType u(2);
    u[0] = capd::IVector(dim); u[0][0] = 1.0;
    u[1] = capd::IVector(dim); u[1][0] = 2.0;

    capd::ddes::DiscreteDelaysFunctionalMap<MockMap, MockSolutionCurve>::ValueStorageType coeffs(2); // Order 1 (coeffs[0], coeffs[1])

    funcMap.computeDDECoefficients(0.0, u, coeffs);

    // Check results
    // coeffs[0] = u[0] = 1.0
    // coeffs[1] = f(u) / 1
    // f(u) = u[0] + u[1] = 1 + 2 = 3.
    // So coeffs[1] should be 3.

    BOOST_CHECK_CLOSE(coeffs[0][0].mid(), 1.0, 1e-10);
    BOOST_CHECK_CLOSE(coeffs[1][0].mid(), 3.0, 1e-10);
}

BOOST_AUTO_TEST_SUITE_END()
