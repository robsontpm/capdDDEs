#define BOOST_TEST_MODULE BasicDiscreteDelaysFunctionalMapBUGGetMaxDelay
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <capd/ddes/BasicDiscreteDelaysFunctionalMap.hpp>
#include <capd/capdlib.h>
#include <vector>

// Mock types (Same as in main test file)
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
            } else {
               coeffs.resize(new_size, VectorType(dimension()));
            }
        }
        return coeffs[k];
    }
};

// Mock Solution Curve
struct MockSolutionCurve {    
    using RealType = ScalarType;
    using VectorType = capd::DVector;
    using MatrixType = capd::DMatrix;
    using TimePointType = double;
    using JetType = MockJet;
    using size_type = std::size_t;
    using DataType = VectorType;

    int dim;

    explicit MockSolutionCurve(int d) : dim(d) {}

    int dimension() const { return dim; }

    VectorType eval(TimePointType t) const {
        VectorType v(dim);
        for(int i=0; i<dim; ++i) v[i] = t;
        return v;
    }

    JetType jet(TimePointType t) const {
        std::vector<VectorType> c(3);
        c[0] = eval(t);
        c[1] = VectorType(dim); for(int i=0; i<dim; ++i) c[1][i] = 1.0;
        c[2] = VectorType(dim);
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
        for(int i=0; i<imgDim; ++i) {
            out[i] = T(0.0);
            for(int j=0; j<in.dimension(); ++j) {
                out[i] = out[i] + in[j];
            }
        }
    }
};

// Test Suite
BOOST_AUTO_TEST_SUITE(BasicDiscreteDelaysFunctionalMapBUGTests)

BOOST_AUTO_TEST_CASE(GetMaxDelayTest) {
    BOOST_WARN_MESSAGE(false, "GetMaxDelayTest is disabled due to compilation error in BasicDiscreteDelaysFunctionalMap::getMaxDelay (bug comparing iterator with value). Uncomment #define BUG_REPRODUCTION to reproduce.");

#ifdef BUG_REPRODUCTION
    int dim = 1;
    MockMap map(dim, dim * 3); // 2 delays
    std::vector<double> delayVec = {1.0, 2.5};

    capd::ddes::BasicDiscreteDelaysFunctionalMap<MockMap, MockSolutionCurve> funcMap(map, delayVec);

    // This line causes compilation error:
    BOOST_CHECK_EQUAL(funcMap.getMaxDelay(), 2.5);
#endif
}

BOOST_AUTO_TEST_SUITE_END()
